function stepA8_behavior_fit_eval(run_id, varargin)
% stepA8_behavior_fit_eval  Fit behavioral models on TRAIN and evaluate on VALID.
%
% Fits behavioral model variants using A7 door-level datasets.
%
% Model 0 (direct trust-as-probability):
%   p_follow = clamp(tau_decision, 0, 1)
%
% Model 1 (baseline threshold):
%   p_follow = sigmoid( k * (tau_decision - self_confidence) )
%
% Model 2 (offset + lapse):
%   z = k * tau_decision + beta * (self_confidence - 0.5)
%   p*_follow = sigmoid(z)
%   p_follow  = (1-eps)*p*_follow + eps*0.5
%
% IMPORTANT (per request): "override" (NOT follow) is treated as the POSITIVE class.
% We therefore also compute probabilities for override:
%   p_override = 1 - p_follow
%
% Outputs are written to:
%   derived/analysis_runs/<run_id>/stepA8_behavior_fit_eval/
%
% Artifacts:
%   - fit_params.mat (MLE estimates, profile likelihood grids, bootstrap samples)
%   - valid_metrics.csv / .mat (pooled + participant-level)
%   - valid_metrics_delta.csv (Δ metrics vs baseline(s))
%   - calibration_bins_valid.csv
%   - A8_plot_data.mat (single source of truth for A9/A6-like reporting)
%   - figures/*.png
%
% Additions (per request):
%   (1) Explicit ΔNLL / ΔBrier / ΔAcc / ΔAUC relative to baseline_global (VALID).
%   (2) Bootstrap distribution plots (k / beta / eps) (TRAIN bootstrap).
%   (3) Additional diagnostic plot for Model 2: p vs linear predictor z.
%   (4) Add Model 0: p = tau_decision as the simplest behavioral mapping.
%   (5) Add random guesser baseline: p_follow=0.5 (VALID) and store Δ vs random_guesser.
%   (6) Add override-focused metrics (override is POSITIVE):
%         - Balanced accuracy (mean recall of follow and override)
%         - Recall_override (TPR for overrides)
%         - Precision_override
%         - F1_override
%         - PR-AUC_override (optional; computed if DoPRAUC=true and both classes present)
%         - NLL_mean_follow / NLL_mean_override (class-conditional NLL)
%         - Brier_follow / Brier_override (class-conditional Brier)
%         - OverrideAcc@0.5 (accuracy where "positive"=override using p_override>=0.5)
%
% Notes:
%   - "baseline_global" uses TRAIN follow rate (mean(y_train)) as constant p_follow on VALID.
%   - "random_guesser_50_50" uses p_follow=0.5 constant predictor on VALID (uninformed).
%   - Δ tables include delta columns vs BOTH baselines.
%
% CHANGE (Jan 2026):
%   - Dropped "baseline_participant_trainrate" baseline entirely, because VALID
%     participants are disjoint from TRAIN in this pipeline, making that baseline
%     degenerate (falls back to baseline_global for all VALID rows).

    if nargin < 1 || isempty(run_id)
        error("stepA8_behavior_fit_eval: run_id is required.");
    end
    run_id = string(run_id);

    p = inputParser;
    p.addParameter("OutDir", "", @(s) isstring(s) || ischar(s));
    p.addParameter("Overwrite", false, @(x) islogical(x) && isscalar(x));
    p.addParameter("NBins", 10, @(x) isnumeric(x) && isscalar(x) && x>=2);
    p.addParameter("BootstrapN", 2000, @(x) isnumeric(x) && isscalar(x) && x>=0);
    p.addParameter("DoAUC", true, @(x) islogical(x) && isscalar(x));
    p.addParameter("DoPRAUC", true, @(x) islogical(x) && isscalar(x)); % PR-AUC for override
    p.addParameter("ProfileGridK", linspace(0, 40, 401), @(x) isnumeric(x) && isvector(x));
    p.addParameter("ProfileGridBeta", linspace(-10, 10, 401), @(x) isnumeric(x) && isvector(x));
    p.addParameter("ProfileGridEps", linspace(0, 0.5, 251), @(x) isnumeric(x) && isvector(x));
    p.addParameter("RandomSeed", 1, @(x) isnumeric(x) && isscalar(x));
    p.parse(varargin{:});
    args = p.Results;

    rng(args.RandomSeed);

    % -------------------------
    % Locate A7 datasets
    % -------------------------
    a7Dir = fullfile("derived","analysis_runs",run_id,"stepA7_behavior_dataset");
    trainMat = fullfile(a7Dir, "behavior_dataset_train.mat");
    validMat = fullfile(a7Dir, "behavior_dataset_valid.mat");
    must_exist_file(trainMat, "A7 TRAIN dataset");
    must_exist_file(validMat, "A7 VALID dataset");

    S_tr = load(trainMat, "T");
    S_va = load(validMat, "T");
    if ~isfield(S_tr,"T") || ~istable(S_tr.T), error("[A8] TRAIN mat missing table T."); end
    if ~isfield(S_va,"T") || ~istable(S_va.T), error("[A8] VALID mat missing table T."); end
    Ttr = S_tr.T;
    Tva = S_va.T;

    % Filter to rows with usable labels
    Ttr = Ttr(Ttr.is_valid_label==1, :);
    Tva = Tva(Tva.is_valid_label==1, :);

    % Minimal required columns
    reqCols = ["participant_id","tau_decision","self_confidence","sc_centered","margin_treshold","followed","block_index","door_index"];
    assert(all(ismember(reqCols, string(Ttr.Properties.VariableNames))), "[A8] TRAIN missing required columns.");
    assert(all(ismember(reqCols, string(Tva.Properties.VariableNames))), "[A8] VALID missing required columns.");

    % Remove NaNs in predictors
    Ttr = Ttr(isfinite(Ttr.tau_decision) & isfinite(Ttr.self_confidence) & isfinite(Ttr.sc_centered) & isfinite(Ttr.margin_treshold), :);
    Tva = Tva(isfinite(Tva.tau_decision) & isfinite(Tva.self_confidence) & isfinite(Tva.sc_centered) & isfinite(Tva.margin_treshold), :);

    % -------------------------
    % Output directory
    % -------------------------
    outDir = string(args.OutDir);
    if strlength(outDir)==0
        outDir = fullfile("derived","analysis_runs",run_id,"stepA8_behavior_fit_eval");
    end
    ensure_dir(outDir);
    figDir = fullfile(outDir, "figures");
    ensure_dir(figDir);

    fitMat     = fullfile(outDir, "fit_params.mat");
    plotMat    = fullfile(outDir, "A8_plot_data.mat");
    metricsCsv = fullfile(outDir, "valid_metrics.csv");
    metricsDeltaCsv = fullfile(outDir, "valid_metrics_delta.csv");
    calibCsv   = fullfile(outDir, "calibration_bins_valid.csv");
    metaMat    = fullfile(outDir, "meta.mat");
    metaJson   = fullfile(outDir, "meta.json");

    if ~args.Overwrite
        if isfile(fitMat) || isfile(metricsCsv) || isfile(plotMat)
            error("[A8] Outputs exist. Set Overwrite=true to replace. (%s)", outDir);
        end
    end

    % -------------------------
    % Prepare arrays
    % -------------------------
    ytr = double(Ttr.followed(:)); % 1=follow, 0=override
    yva = double(Tva.followed(:));

    % Override labels (positive class)
    ytr_ov = 1 - ytr; % 1=override
    yva_ov = 1 - yva;

    m_tr = double(Ttr.margin_treshold(:)); % tau - sc
    m_va = double(Tva.margin_treshold(:));

    tau_tr = double(Ttr.tau_decision(:));
    tau_va = double(Tva.tau_decision(:));

    scC_tr = double(Ttr.sc_centered(:));   % sc - 0.5
    scC_va = double(Tva.sc_centered(:));

    pid_tr = string(Ttr.participant_id); %#ok<NASGU> % kept for bootstrap resampling
    pid_va = string(Tva.participant_id);

    % -------------------------
    % Model 0: Direct p_follow = tau (no fitting)
    % -------------------------
    p0_va_follow = clamp01(tau_va);
    p0_va_ov     = 1 - p0_va_follow;

    % -------------------------
    % Fit Model 1: k (MLE)
    % -------------------------
    f1 = @(z) nll_model1(exp(z), m_tr, ytr);
    z0_1 = log(10);
    zhat1 = fminsearch(f1, z0_1, optimset('Display','off'));
    k1_hat = exp(zhat1);

    kGrid = args.ProfileGridK(:);
    nll1_grid = NaN(numel(kGrid),1);
    for i = 1:numel(kGrid)
        nll1_grid(i) = nll_model1(kGrid(i), m_tr, ytr);
    end

    % -------------------------
    % Fit Model 2: (k, beta, eps) (MLE)
    % -------------------------
    epsMax = 0.5;
    f2 = @(z) nll_model2(exp(z(1)), z(2), epsMax*sigmoid(z(3)), tau_tr, scC_tr, ytr);

    z0_2 = [log(10); 0; logit(0.05/epsMax)];
    zhat2 = fminsearch(f2, z0_2, optimset('Display','off'));
    k2_hat    = exp(zhat2(1));
    beta_hat  = zhat2(2);
    eps_hat   = epsMax*sigmoid(zhat2(3));

    betaGrid = args.ProfileGridBeta(:);
    epsGrid  = args.ProfileGridEps(:);

    nll2_kgrid = NaN(numel(kGrid),1);
    for i = 1:numel(kGrid)
        nll2_kgrid(i) = nll_model2(kGrid(i), beta_hat, eps_hat, tau_tr, scC_tr, ytr);
    end

    nll2_betagrid = NaN(numel(betaGrid),1);
    for i = 1:numel(betaGrid)
        nll2_betagrid(i) = nll_model2(k2_hat, betaGrid(i), eps_hat, tau_tr, scC_tr, ytr);
    end

    nll2_epsgrid = NaN(numel(epsGrid),1);
    for i = 1:numel(epsGrid)
        nll2_epsgrid(i) = nll_model2(k2_hat, beta_hat, epsGrid(i), tau_tr, scC_tr, ytr);
    end

    % -------------------------
    % Bootstrap uncertainty (participant-level)
    % -------------------------
    B = args.BootstrapN;
    boot = struct();
    boot.B = B;

    if B > 0
        uniqP = unique(string(Ttr.participant_id));
        nP = numel(uniqP);

        boot.k1 = NaN(B,1);
        boot.k2 = NaN(B,1);
        boot.beta = NaN(B,1);
        boot.eps = NaN(B,1);

        for b = 1:B
            sampP = uniqP(randi(nP, nP, 1)); % resample participants w/ replacement

            % Build indices with replication (correct bootstrap weighting)
            idx = [];
            for sp = 1:numel(sampP)
                idx_sp = find(string(Ttr.participant_id) == sampP(sp));
                idx = [idx; idx_sp]; %#ok<AGROW>
            end

            yb = ytr(idx);
            mb = m_tr(idx);
            taub = tau_tr(idx);
            scCb = scC_tr(idx);

            % Fit Model 1
            f1b = @(z) nll_model1(exp(z), mb, yb);
            zb = fminsearch(f1b, z0_1, optimset('Display','off'));
            boot.k1(b) = exp(zb);

            % Fit Model 2
            f2b = @(z) nll_model2(exp(z(1)), z(2), epsMax*sigmoid(z(3)), taub, scCb, yb);
            zb2 = fminsearch(f2b, z0_2, optimset('Display','off'));
            boot.k2(b)   = exp(zb2(1));
            boot.beta(b) = zb2(2);
            boot.eps(b)  = epsMax*sigmoid(zb2(3));
        end

        boot.ci95_k1 = prctile(boot.k1, [2.5 97.5]);
        boot.ci95_k2   = prctile(boot.k2, [2.5 97.5]);
        boot.ci95_beta = prctile(boot.beta, [2.5 97.5]);
        boot.ci95_eps  = prctile(boot.eps, [2.5 97.5]);
    end

    % -------------------------
    % VALID predictions
    % -------------------------
    % Model 1
    p1_va_follow = sigmoid(k1_hat .* m_va);
    p1_va_ov     = 1 - p1_va_follow;

    % Model 2
    z2_va = k2_hat.*tau_va + beta_hat.*scC_va;
    p2_va_follow = (1-eps_hat).*sigmoid(z2_va) + eps_hat.*0.5;
    p2_va_follow = clamp01(p2_va_follow);
    p2_va_ov     = 1 - p2_va_follow;

    % Baseline predictors
    p_base_global_follow = mean(ytr); % train follow-rate
    p_base_global_va = p_base_global_follow * ones(size(yva));
    p_base_global_ov = 1 - p_base_global_va;

    % Random guesser (uninformed)
    p_rand_follow_va = 0.5 * ones(size(yva));
    p_rand_ov_va     = 0.5 * ones(size(yva));

    % -------------------------
    % VALID pooled metrics (override-focused metrics included)
    % -------------------------
    rows = {
        "random_guesser_50_50", metrics_all_with_override(yva, yva_ov, p_rand_follow_va, p_rand_ov_va, args.DoAUC, args.DoPRAUC);
        "baseline_global", metrics_all_with_override(yva, yva_ov, p_base_global_va, p_base_global_ov, args.DoAUC, args.DoPRAUC);
        "model0_trust_as_probability", metrics_all_with_override(yva, yva_ov, p0_va_follow, p0_va_ov, args.DoAUC, args.DoPRAUC);
        "model1_threshold", metrics_all_with_override(yva, yva_ov, p1_va_follow, p1_va_ov, args.DoAUC, args.DoPRAUC);
        "model2_offset_lapse", metrics_all_with_override(yva, yva_ov, p2_va_follow, p2_va_ov, args.DoAUC, args.DoPRAUC);
    };

    validOverall = cell2table(rows, 'VariableNames', {'method','metrics'});
    validOverall = unpack_metrics(validOverall);

    % Participant-level metrics (VALID)
    validByP = metrics_by_participant(pid_va, yva, yva_ov, ...
        p0_va_follow, p0_va_ov, ...
        p1_va_follow, p1_va_ov, ...
        p2_va_follow, p2_va_ov, ...
        p_base_global_va, p_base_global_ov, ...
        p_rand_follow_va, p_rand_ov_va, ...
        args.DoAUC, args.DoPRAUC);

    % Calibration bins (VALID) for Model 0/1/2 (follow-prob calibration)
    calib0 = calibration_bins(yva, p0_va_follow, args.NBins);
    calib0.method = repmat("model0_trust_as_probability", height(calib0), 1);

    calib1 = calibration_bins(yva, p1_va_follow, args.NBins);
    calib1.method = repmat("model1_threshold", height(calib1), 1);

    calib2 = calibration_bins(yva, p2_va_follow, args.NBins);
    calib2.method = repmat("model2_offset_lapse", height(calib2), 1);

    calib = [calib0; calib1; calib2];

    writetable(validOverall, metricsCsv);
    writetable(calib, calibCsv);

    % -------------------------
    % Δ metrics vs baselines (VALID)
    % -------------------------
    validDelta = make_delta_table_multi(validOverall, ["baseline_global","random_guesser_50_50"]);
    writetable(validDelta, metricsDeltaCsv);

    % -------------------------
    % Save MAT bundles + meta
    % -------------------------
    fit = struct();

    fit.model0 = struct();
    fit.model0.description = "p_follow = tau_decision (clamped to [0,1])";

    fit.model1.k_hat = k1_hat;
    fit.model1.k_grid = kGrid;
    fit.model1.nll_grid = nll1_grid;

    fit.model2.k_hat = k2_hat;
    fit.model2.beta_hat = beta_hat;
    fit.model2.eps_hat = eps_hat;
    fit.model2.k_grid = kGrid;
    fit.model2.nll_kgrid = nll2_kgrid;
    fit.model2.beta_grid = betaGrid;
    fit.model2.nll_betagrid = nll2_betagrid;
    fit.model2.eps_grid = epsGrid;
    fit.model2.nll_epsgrid = nll2_epsgrid;

    fit.bootstrap = boot;

    save(fitMat, "fit", "validOverall", "validDelta", "validByP", "calib", "-v7.3");

    A8plot = struct();
    A8plot.run_id = char(run_id);
    A8plot.fit = fit;
    A8plot.validOverall = validOverall;
    A8plot.validDelta = validDelta;
    A8plot.validByParticipant = validByP;
    A8plot.calibration = calib;

    A8plot.pred_valid = struct( ...
        "y_follow",yva, ...
        "y_override",yva_ov, ...
        "p_random_follow",p_rand_follow_va, ...
        "p_random_override",p_rand_ov_va, ...
        "p_model0_follow",p0_va_follow, ...
        "p_model0_override",p0_va_ov, ...
        "p_model1_follow",p1_va_follow, ...
        "p_model1_override",p1_va_ov, ...
        "p_model2_follow",p2_va_follow, ...
        "p_model2_override",p2_va_ov, ...
        "z_model2", z2_va, ...
        "p_base_global_follow",p_base_global_va, ...
        "p_base_global_override",p_base_global_ov, ...
        "pid",pid_va, ...
        "block", double(Tva.block_index(:)), ...
        "door", double(Tva.door_index(:)), ...
        "tau", tau_va, ...
        "sc_centered", scC_va, ...
        "margin_treshold", m_va);

    save(plotMat, "A8plot", "-v7.3");

    meta = struct();
    meta.run_id = char(run_id);
    meta.created = char(datetime('now','Format','yyyy-MM-dd HH:mm:ss'));
    meta.train_rows = height(Ttr);
    meta.valid_rows = height(Tva);
    meta.bootstrapN = B;
    meta.nbins = args.NBins;
    meta.random_seed = args.RandomSeed;
    meta.override_positive = true;
    meta.do_auc = args.DoAUC;
    meta.do_prauc = args.DoPRAUC;
    meta.dropped_participant_trainrate_baseline = true;
    save(metaMat, "meta");
    save_json(metaJson, meta);

    % -------------------------
    % Figures
    % -------------------------
    % Profiles (TRAIN)
    make_fig_profile_k(fullfile(figDir, "train_profile_model1_k.png"), kGrid, nll1_grid, k1_hat, "Model 1: profile NLL(k)");
    make_fig_profile_k(fullfile(figDir, "train_profile_model2_k.png"), kGrid, nll2_kgrid, k2_hat, "Model 2: profile NLL(k) (others fixed)");
    make_fig_profile_1d(fullfile(figDir, "train_profile_model2_beta.png"), betaGrid, nll2_betagrid, beta_hat, "Model 2: profile NLL(beta)");
    make_fig_profile_1d(fullfile(figDir, "train_profile_model2_eps.png"), epsGrid, nll2_epsgrid, eps_hat, "Model 2: profile NLL(eps)");

    % Calibration (VALID) - for follow probability
    make_fig_calibration(fullfile(figDir, "valid_calibration_model0.png"), calib0, "VALID calibration: Model 0 (p_follow=tau)");
    make_fig_calibration(fullfile(figDir, "valid_calibration_model1.png"), calib1, "VALID calibration: Model 1");
    make_fig_calibration(fullfile(figDir, "valid_calibration_model2.png"), calib2, "VALID calibration: Model 2");

    % p vs x (VALID)
    make_fig_prob_vs_x(fullfile(figDir, "valid_prob_vs_tau_model0.png"), tau_va, p0_va_follow, yva, ...
        "VALID: p_follow vs tau (Model 0)", "tau\_decision");

    make_fig_prob_vs_x(fullfile(figDir, "valid_prob_vs_margin_model1.png"), m_va, p1_va_follow, yva, ...
        "VALID: p_follow vs (tau-sc) (Model 1)", "margin (tau-sc)");

    make_fig_prob_vs_x(fullfile(figDir, "valid_prob_vs_margin_model2.png"), m_va, p2_va_follow, yva, ...
        "VALID: p_follow vs (tau-sc) (Model 2, plotted vs same margin)", "margin (tau-sc)");

    % Model 2: p vs z
    make_fig_prob_vs_x(fullfile(figDir, "valid_prob_vs_z_model2.png"), z2_va, p2_va_follow, yva, ...
        "VALID: p_follow vs z (Model 2)", "z = k*tau + beta*(sc-0.5)");

    % PR curve for override (positive class) for the three models
    if args.DoPRAUC
        make_fig_pr_curve(fullfile(figDir, "valid_prcurve_override_models012.png"), ...
            yva_ov, ...
            {p0_va_ov, p1_va_ov, p2_va_ov}, ...
            ["Model0","Model1","Model2"]);
    end

    % NLL-by-class plot for the three models (follow vs override)
    make_fig_class_nll(fullfile(figDir, "valid_nll_by_class_models012.png"), validOverall);

    % Bootstrap distributions (TRAIN)
    if B > 0
        make_fig_bootstrap_hists(fullfile(figDir, "train_bootstrap_params_model1.png"), ...
            boot.k1, "Bootstrap: Model 1 (k)", "k");
        make_fig_bootstrap_hists(fullfile(figDir, "train_bootstrap_params_model2.png"), ...
            [boot.k2, boot.beta, boot.eps], "Bootstrap: Model 2 (k, beta, eps)", ["k","beta","eps"]);
    end

    fprintf("[Step A8] Done.\n");
    fprintf("  Model0: p_follow=tau (no fit)\n");
    fprintf("  Model1: k=%.4g\n", k1_hat);
    fprintf("  Model2: k=%.4g, beta=%.4g, eps=%.4g\n", k2_hat, beta_hat, eps_hat);
    fprintf("  Output dir: %s\n", outDir);

    % Helpful terminal summary for ΔNLL_mean vs baselines + override recall quick glance
    try
        show_delta_summary_multi(validDelta, "baseline_global", "random_guesser_50_50");
        show_override_quick_summary(validOverall);
    catch
        % no-op
    end
end

% ======================================================================
% Helpers: NLL
% ======================================================================

function nll = nll_model1(k, margin, y_follow)
    if ~isfinite(k) || k < 0, nll = Inf; return; end
    p = sigmoid(k .* margin);
    nll = bernoulli_nll(y_follow, p);
end

function nll = nll_model2(k, beta, eps, tau, scC, y_follow)
    if ~isfinite(k) || k < 0 || ~isfinite(beta) || ~isfinite(eps) || eps < 0 || eps > 1
        nll = Inf; return;
    end
    z = k .* tau + beta .* scC;
    pstar = sigmoid(z);
    p = (1-eps).*pstar + eps.*0.5;
    nll = bernoulli_nll(y_follow, p);
end

function nll = bernoulli_nll(y, p)
    p = min(max(p, 1e-12), 1-1e-12);
    nll = -sum(y .* log(p) + (1-y).*log(1-p));
end

function p = clamp01(x)
    p = min(max(x, 0), 1);
end

% ======================================================================
% Helpers: metrics (pooled)
% ======================================================================

function M = metrics_all_with_override(y_follow, y_override, p_follow, p_override, doAUC, doPRAUC)
    % y_follow: 1=follow, 0=override
    % y_override: 1=override, 0=follow  (positive class)
    %
    % p_follow: predicted probability of follow
    % p_override: predicted probability of override (should be 1-p_follow)

    p_follow   = min(max(p_follow, 1e-12), 1-1e-12);
    p_override = min(max(p_override, 1e-12), 1-1e-12);

    n = numel(y_follow);

    % Proper scoring rules (for follow label as originally defined)
    nll_total = -sum(y_follow .* log(p_follow) + (1-y_follow).*log(1-p_follow));
    brier = mean((p_follow - y_follow).^2);

    % Old "accuracy" (threshold 0.5 on follow prob)
    acc_follow = mean((p_follow >= 0.5) == (y_follow==1));

    % Override classification at 0.5 threshold (positive=override)
    pred_ov = (p_override >= 0.5);
    [prec_ov, rec_ov, f1_ov, bal_acc, acc_ovpos] = override_metrics_from_predictions(y_override, pred_ov);

    % Class-conditional scoring
    [nll_mean_follow, nll_mean_override] = nll_by_class(y_follow, p_follow);
    [brier_follow, brier_override] = brier_by_class(y_follow, p_follow);

    out = struct();
    out.N = n;

    out.NLL_total = nll_total;
    out.NLL_mean  = nll_total / max(1,n);
    out.Brier     = brier;
    out.Acc       = acc_follow;

    % Original ROC-AUC (for follow as positive) if requested
    if doAUC
        out.AUC = auc_roc(y_follow, p_follow);
    else
        out.AUC = NaN;
    end

    % Override-focused metrics (override is positive)
    out.OverrideRate = mean(y_override==1);
    out.OverridePrecision = prec_ov;
    out.OverrideRecall = rec_ov;
    out.OverrideF1 = f1_ov;
    out.BalancedAcc = bal_acc;
    out.OverrideAcc_at_0_5 = acc_ovpos;

    if doPRAUC
        out.OverridePRAUC = pr_auc(y_override, p_override);
    else
        out.OverridePRAUC = NaN;
    end

    out.NLL_mean_follow = nll_mean_follow;
    out.NLL_mean_override = nll_mean_override;
    out.Brier_follow = brier_follow;
    out.Brier_override = brier_override;

    M = out;
end

function [prec_ov, rec_ov, f1_ov, bal_acc, acc_ovpos] = override_metrics_from_predictions(y_override, pred_override)
    % y_override: logical/0-1
    y = (y_override==1);
    p = (pred_override==1);

    TP = sum(p & y);
    FP = sum(p & ~y);
    FN = sum(~p & y);
    TN = sum(~p & ~y);

    % Precision/Recall for override (positive class)
    prec_ov = TP / max(1, (TP+FP));
    rec_ov  = TP / max(1, (TP+FN));
    f1_ov   = (2*prec_ov*rec_ov) / max(1e-12, (prec_ov+rec_ov));

    % Recall for follow is TNR here (since follow is negative class)
    rec_follow = TN / max(1, (TN+FP));

    bal_acc = 0.5*(rec_ov + rec_follow);

    % "Accuracy" when treating override as the positive class
    acc_ovpos = (TP + TN) / max(1, (TP+FP+FN+TN));
end

function [nll_f, nll_o] = nll_by_class(y_follow, p_follow)
    % Mean NLL over follow trials and over override trials.
    p_follow = min(max(p_follow, 1e-12), 1-1e-12);
    maskF = (y_follow==1);
    maskO = (y_follow==0);

    if any(maskF)
        nll_f = mean(-log(p_follow(maskF)));
    else
        nll_f = NaN;
    end
    if any(maskO)
        nll_o = mean(-log(1 - p_follow(maskO)));
    else
        nll_o = NaN;
    end
end

function [b_f, b_o] = brier_by_class(y_follow, p_follow)
    maskF = (y_follow==1);
    maskO = (y_follow==0);
    if any(maskF)
        b_f = mean((p_follow(maskF) - 1).^2);
    else
        b_f = NaN;
    end
    if any(maskO)
        b_o = mean((p_follow(maskO) - 0).^2);
    else
        b_o = NaN;
    end
end

% ======================================================================
% Helpers: unpack metrics
% ======================================================================

function T = unpack_metrics(Tin)
    methods = Tin.method;
    S = Tin.metrics;

    if iscell(S)
        s1 = S{1};
        getField = @(r,f) S{r}.(f);
    elseif isstruct(S)
        s1 = S(1);
        getField = @(r,f) S(r).(f);
    else
        error("unpack_metrics: Unexpected type for Tin.metrics: %s", class(S));
    end

    fn = fieldnames(s1);

    T = table();
    T.method = methods;

    for i = 1:numel(fn)
        f = fn{i};
        vals = NaN(height(Tin),1);
        for r = 1:height(Tin)
            vals(r) = getField(r,f);
        end
        T.(f) = vals;
    end
end

% ======================================================================
% Helpers: participant-level metrics (NO participant-trainrate baseline)
% ======================================================================

function Tp = metrics_by_participant(pid, y_follow, y_override, ...
    p0_follow, p0_ov, ...
    p1_follow, p1_ov, ...
    p2_follow, p2_ov, ...
    pb_follow, pb_ov, ...
    pr_follow, pr_ov, ...
    doAUC, doPRAUC)

    uniqP = unique(pid);
    Tp = table();
    Tp.participant_id = uniqP;
    Tp.N = zeros(numel(uniqP),1);

    cols = [ ...
        "OverrideRate", ...
        "NLL_mean_random","Brier_random","Acc_random","AUC_random","OverrideRecall_random","OverridePrecision_random","OverrideF1_random","BalancedAcc_random","OverridePRAUC_random","NLL_mean_override_random", ...
        "NLL_mean_base_global","Brier_base_global","Acc_base_global","AUC_base_global","OverrideRecall_base_global","OverridePrecision_base_global","OverrideF1_base_global","BalancedAcc_base_global","OverridePRAUC_base_global","NLL_mean_override_base_global", ...
        "NLL_mean_model0","Brier_model0","Acc_model0","AUC_model0","OverrideRecall_model0","OverridePrecision_model0","OverrideF1_model0","BalancedAcc_model0","OverridePRAUC_model0","NLL_mean_override_model0", ...
        "NLL_mean_model1","Brier_model1","Acc_model1","AUC_model1","OverrideRecall_model1","OverridePrecision_model1","OverrideF1_model1","BalancedAcc_model1","OverridePRAUC_model1","NLL_mean_override_model1", ...
        "NLL_mean_model2","Brier_model2","Acc_model2","AUC_model2","OverrideRecall_model2","OverridePrecision_model2","OverrideF1_model2","BalancedAcc_model2","OverridePRAUC_model2","NLL_mean_override_model2" ...
    ];

    for c = 1:numel(cols)
        Tp.(cols(c)) = NaN(numel(uniqP),1);
    end

    for i = 1:numel(uniqP)
        mask = (pid==uniqP(i));
        Tp.N(i) = sum(mask);

        Tp.OverrideRate(i) = mean(y_override(mask)==1);

        mr  = metrics_all_with_override(y_follow(mask), y_override(mask), pr_follow(mask), pr_ov(mask), doAUC, doPRAUC);
        mbg = metrics_all_with_override(y_follow(mask), y_override(mask), pb_follow(mask), pb_ov(mask), doAUC, doPRAUC);

        mm0 = metrics_all_with_override(y_follow(mask), y_override(mask), p0_follow(mask), p0_ov(mask), doAUC, doPRAUC);
        mm1 = metrics_all_with_override(y_follow(mask), y_override(mask), p1_follow(mask), p1_ov(mask), doAUC, doPRAUC);
        mm2 = metrics_all_with_override(y_follow(mask), y_override(mask), p2_follow(mask), p2_ov(mask), doAUC, doPRAUC);

        % Random
        Tp.NLL_mean_random(i) = mr.NLL_mean;
        Tp.Brier_random(i)    = mr.Brier;
        Tp.Acc_random(i)      = mr.Acc;
        Tp.AUC_random(i)      = mr.AUC;
        Tp.OverrideRecall_random(i)    = mr.OverrideRecall;
        Tp.OverridePrecision_random(i) = mr.OverridePrecision;
        Tp.OverrideF1_random(i)        = mr.OverrideF1;
        Tp.BalancedAcc_random(i)       = mr.BalancedAcc;
        Tp.OverridePRAUC_random(i)      = mr.OverridePRAUC;
        Tp.NLL_mean_override_random(i)  = mr.NLL_mean_override;

        % Base global
        Tp.NLL_mean_base_global(i) = mbg.NLL_mean;
        Tp.Brier_base_global(i)    = mbg.Brier;
        Tp.Acc_base_global(i)      = mbg.Acc;
        Tp.AUC_base_global(i)      = mbg.AUC;
        Tp.OverrideRecall_base_global(i)    = mbg.OverrideRecall;
        Tp.OverridePrecision_base_global(i) = mbg.OverridePrecision;
        Tp.OverrideF1_base_global(i)        = mbg.OverrideF1;
        Tp.BalancedAcc_base_global(i)       = mbg.BalancedAcc;
        Tp.OverridePRAUC_base_global(i)      = mbg.OverridePRAUC;
        Tp.NLL_mean_override_base_global(i)  = mbg.NLL_mean_override;

        % Model 0
        Tp.NLL_mean_model0(i) = mm0.NLL_mean;
        Tp.Brier_model0(i)    = mm0.Brier;
        Tp.Acc_model0(i)      = mm0.Acc;
        Tp.AUC_model0(i)      = mm0.AUC;
        Tp.OverrideRecall_model0(i)    = mm0.OverrideRecall;
        Tp.OverridePrecision_model0(i) = mm0.OverridePrecision;
        Tp.OverrideF1_model0(i)        = mm0.OverrideF1;
        Tp.BalancedAcc_model0(i)       = mm0.BalancedAcc;
        Tp.OverridePRAUC_model0(i)      = mm0.OverridePRAUC;
        Tp.NLL_mean_override_model0(i)  = mm0.NLL_mean_override;

        % Model 1
        Tp.NLL_mean_model1(i) = mm1.NLL_mean;
        Tp.Brier_model1(i)    = mm1.Brier;
        Tp.Acc_model1(i)      = mm1.Acc;
        Tp.AUC_model1(i)      = mm1.AUC;
        Tp.OverrideRecall_model1(i)    = mm1.OverrideRecall;
        Tp.OverridePrecision_model1(i) = mm1.OverridePrecision;
        Tp.OverrideF1_model1(i)        = mm1.OverrideF1;
        Tp.BalancedAcc_model1(i)       = mm1.BalancedAcc;
        Tp.OverridePRAUC_model1(i)      = mm1.OverridePRAUC;
        Tp.NLL_mean_override_model1(i)  = mm1.NLL_mean_override;

        % Model 2
        Tp.NLL_mean_model2(i) = mm2.NLL_mean;
        Tp.Brier_model2(i)    = mm2.Brier;
        Tp.Acc_model2(i)      = mm2.Acc;
        Tp.AUC_model2(i)      = mm2.AUC;
        Tp.OverrideRecall_model2(i)    = mm2.OverrideRecall;
        Tp.OverridePrecision_model2(i) = mm2.OverridePrecision;
        Tp.OverrideF1_model2(i)        = mm2.OverrideF1;
        Tp.BalancedAcc_model2(i)       = mm2.BalancedAcc;
        Tp.OverridePRAUC_model2(i)      = mm2.OverridePRAUC;
        Tp.NLL_mean_override_model2(i)  = mm2.NLL_mean_override;
    end
end

% ======================================================================
% Helpers: calibration
% ======================================================================

function C = calibration_bins(y, p, nbins)
    p = min(max(p, 0), 1);

    edges = linspace(0,1,nbins+1);
    bin = discretize(p, edges);
    bin(isnan(bin)) = nbins;

    C = table();
    C.bin = (1:nbins)';
    C.p_lo = edges(1:end-1)';
    C.p_hi = edges(2:end)';

    C.n = zeros(nbins,1);
    C.p_mean = NaN(nbins,1);
    C.y_mean = NaN(nbins,1);
    C.ci_lo = NaN(nbins,1);
    C.ci_hi = NaN(nbins,1);

    for b = 1:nbins
        mask = (bin==b);
        C.n(b) = sum(mask);
        if C.n(b) > 0
            C.p_mean(b) = mean(p(mask));
            C.y_mean(b) = mean(y(mask));
            [lo, hi] = wilson_ci(sum(y(mask)), C.n(b), 0.05);
            C.ci_lo(b) = lo;
            C.ci_hi(b) = hi;
        end
    end

    w = C.n / max(1,sum(C.n));
    C.abs_gap = abs(C.p_mean - C.y_mean);
    ece = nansum(w .* C.abs_gap);
    C.ECE = repmat(ece, nbins, 1);
end

function [lo, hi] = wilson_ci(k, n, alpha)
    if n <= 0
        lo = NaN; hi = NaN; return;
    end
    z = norminv(1 - alpha/2);
    phat = k/n;
    denom = 1 + z^2/n;
    center = (phat + z^2/(2*n)) / denom;
    half = (z/denom) * sqrt((phat*(1-phat) + z^2/(4*n))/n);
    lo = max(0, center - half);
    hi = min(1, center + half);
end

% ======================================================================
% Helpers: ROC-AUC and PR-AUC
% ======================================================================

function a = auc_roc(y, p)
    y = y(:);
    p = p(:);
    if numel(unique(y)) < 2
        a = NaN; return;
    end
    r = tiedrank(p);
    n1 = sum(y==1);
    n0 = sum(y==0);
    a = (sum(r(y==1)) - n1*(n1+1)/2) / (n1*n0);
end

function ap = pr_auc(y_pos, p_pos)
    % Average precision (area under precision-recall curve)
    y_pos = y_pos(:) == 1;
    p_pos = p_pos(:);

    if numel(unique(y_pos)) < 2
        ap = NaN; return;
    end

    [~, idx] = sort(p_pos, 'descend');
    y_sorted = y_pos(idx);

    tp = cumsum(y_sorted==1);
    fp = cumsum(y_sorted==0);

    prec = tp ./ max(1, (tp + fp));
    rec  = tp ./ max(1, sum(y_sorted==1));

    drec = [rec(1); diff(rec)];
    ap = sum(prec .* drec);
    ap = min(max(ap, 0), 1);
end

% ======================================================================
% Helpers: sigmoid/logit
% ======================================================================

function s = sigmoid(z)
    s = 1 ./ (1 + exp(-z));
end

function z = logit(p)
    p = min(max(p, 1e-12), 1-1e-12);
    z = log(p./(1-p));
end

% ======================================================================
% Δ metrics table helper (multi-baseline)
% ======================================================================

function Td = make_delta_table_multi(T, baselineMethods)
    if ~ismember("method", string(T.Properties.VariableNames))
        error("make_delta_table_multi: T must include 'method'.");
    end

    Td = T;
    metricNames = ["NLL_total","NLL_mean","Brier","Acc","AUC", ...
                   "OverridePrecision","OverrideRecall","OverrideF1","BalancedAcc","OverridePRAUC", ...
                   "NLL_mean_override","NLL_mean_follow", "Brier_override","Brier_follow", ...
                   "OverrideAcc_at_0_5"];

    baselineMethods = string(baselineMethods(:))';
    for bm = baselineMethods
        baseIdx = find(string(Td.method) == bm, 1);
        if isempty(baseIdx)
            error("make_delta_table_multi: baseline method '%s' not found.", bm);
        end
        base = Td(baseIdx, :);

        for mn = metricNames
            if ismember(mn, string(Td.Properties.VariableNames))
                Td.("d_" + mn + "_vs_" + bm) = Td.(mn) - base.(mn);
            end
        end
    end
end

function show_delta_summary_multi(Td, baselineGlobalName, randomName)
    keyG = "d_NLL_mean_vs_" + string(baselineGlobalName);
    keyR = "d_NLL_mean_vs_" + string(randomName);

    haveG = ismember(keyG, string(Td.Properties.VariableNames));
    haveR = ismember(keyR, string(Td.Properties.VariableNames));

    rows = ["model0_trust_as_probability","model1_threshold","model2_offset_lapse"];
    labels = ["Model 0","Model 1","Model 2"];

    for i = 1:numel(rows)
        r = Td(string(Td.method)==rows(i), :);
        if isempty(r), continue; end
        if haveG
            fprintf("  ΔNLL_mean vs %s (%s): %+0.4f\n", baselineGlobalName, labels(i), r.(keyG));
        end
        if haveR
            fprintf("  ΔNLL_mean vs %s (%s): %+0.4f\n", randomName, labels(i), r.(keyR));
        end
    end
end

function show_override_quick_summary(validOverall)
    fprintf("  Override-focused (positive=override): Recall / Precision / F1 / BalancedAcc\n");
    rows = ["random_guesser_50_50","baseline_global","model0_trust_as_probability","model1_threshold","model2_offset_lapse"];
    for i = 1:numel(rows)
        r = validOverall(string(validOverall.method)==rows(i), :);
        if isempty(r), continue; end
        fprintf("    %-28s: R=%.3f  P=%.3f  F1=%.3f  BalAcc=%.3f\n", ...
            rows(i), r.OverrideRecall, r.OverridePrecision, r.OverrideF1, r.BalancedAcc);
    end
end

% ======================================================================
% Plot helpers
% ======================================================================

function make_fig_profile_k(pathPng, x, nll, xhat, titleStr)
    f = figure('Visible','off');
    plot(x, nll, 'LineWidth', 1.5);
    hold on;
    yl = ylim;
    plot([xhat xhat], yl, '--', 'LineWidth', 1.0);
    xlabel('parameter');
    ylabel('NLL (TRAIN)');
    title(titleStr, 'Interpreter','none');
    grid on;
    saveas(f, pathPng);
    close(f);
end

function make_fig_profile_1d(pathPng, x, nll, xhat, titleStr)
    f = figure('Visible','off');
    plot(x, nll, 'LineWidth', 1.5);
    hold on;
    yl = ylim;
    plot([xhat xhat], yl, '--', 'LineWidth', 1.0);
    xlabel('parameter');
    ylabel('NLL (TRAIN)');
    title(titleStr, 'Interpreter','none');
    grid on;
    saveas(f, pathPng);
    close(f);
end

function make_fig_calibration(pathPng, C, titleStr)
    f = figure('Visible','off');
    plot([0 1],[0 1], '--', 'LineWidth', 1.0);
    hold on;

    mask = C.n > 0 & isfinite(C.p_mean) & isfinite(C.y_mean);
    x = C.p_mean(mask);
    y = C.y_mean(mask);
    lo = C.ci_lo(mask);
    hi = C.ci_hi(mask);
    errLow = y - lo;
    errHigh = hi - y;
    errorbar(x, y, errLow, errHigh, 'o', 'LineWidth', 1.0);

    xlabel('Predicted probability (bin mean)');
    ylabel('Empirical follow rate');
    title(sprintf('%s (ECE=%.3f)', titleStr, C.ECE(1)), 'Interpreter','none');
    grid on;
    xlim([0 1]); ylim([0 1]);
    saveas(f, pathPng);
    close(f);
end

function make_fig_prob_vs_x(pathPng, x, p, y_follow, titleStr, xlabStr)
    f = figure('Visible','off');

    n = numel(x);
    idx = 1:n;
    if n > 800
        idx = idx(randperm(n, 800));
    end

    scatter(x(idx), p(idx), 12, y_follow(idx), 'filled');
    xlabel(xlabStr);
    ylabel('Predicted p(follow)');
    title(titleStr, 'Interpreter','none');
    grid on;
    saveas(f, pathPng);
    close(f);
end

function make_fig_bootstrap_hists(pathPng, X, titleStr, paramNames)
    f = figure('Visible','off');

    if isvector(X)
        x = X(:);
        x = x(isfinite(x));
        histogram(x, 40);
        xlabel(string(paramNames));
        ylabel('count');
        title(titleStr, 'Interpreter','none');
        grid on;
        saveas(f, pathPng);
        close(f);
        return;
    end

    [~,D] = size(X);
    if ischar(paramNames) || isstring(paramNames)
        paramNames = string(paramNames);
    end
    if numel(paramNames) ~= D
        paramNames = "param" + (1:D);
    end

    for d = 1:D
        subplot(D,1,d);
        xd = X(:,d);
        xd = xd(isfinite(xd));
        histogram(xd, 40);
        xlabel(paramNames(d));
        ylabel('count');
        grid on;
        if d == 1
            title(titleStr, 'Interpreter','none');
        end
    end

    saveas(f, pathPng);
    close(f);
end

function make_fig_pr_curve(pathPng, y_pos, ppos_list, label_list)
    f = figure('Visible','off');
    hold on;
    grid on;

    prev = mean(y_pos==1);
    plot([0 1], [prev prev], '--', 'LineWidth', 1.0);

    for i = 1:numel(ppos_list)
        p = ppos_list{i}(:);
        y = (y_pos(:)==1);

        if numel(unique(y)) < 2
            continue;
        end

        [~, idx] = sort(p, 'descend');
        y_sorted = y(idx);

        tp = cumsum(y_sorted==1);
        fp = cumsum(y_sorted==0);

        prec = tp ./ max(1, (tp + fp));
        rec  = tp ./ max(1, sum(y_sorted==1));

        plot(rec, prec, 'LineWidth', 1.5);
    end

    xlabel('Recall (override)');
    ylabel('Precision (override)');
    title('VALID PR curves (positive = override)', 'Interpreter','none');

    lg = ["prevalence"] + string(label_list);
    legend(lg, 'Interpreter','none', 'Location','best');

    xlim([0 1]); ylim([0 1]);
    saveas(f, pathPng);
    close(f);
end

function make_fig_class_nll(pathPng, validOverall)
    f = figure('Visible','off');
    grid on; hold on;

    methods = ["model0_trust_as_probability","model1_threshold","model2_offset_lapse"];
    xs = 1:numel(methods);

    nllF = NaN(size(xs));
    nllO = NaN(size(xs));

    for i = 1:numel(methods)
        r = validOverall(string(validOverall.method)==methods(i), :);
        if isempty(r), continue; end
        nllF(i) = r.NLL_mean_follow;
        nllO(i) = r.NLL_mean_override;
    end

    plot(xs, nllF, 'o-', 'LineWidth', 1.5);
    plot(xs, nllO, 's-', 'LineWidth', 1.5);

    set(gca,'XTick',xs,'XTickLabel',methods);
    xtickangle(20);
    ylabel('Mean NLL (lower is better)');
    title('VALID class-conditional NLL (follow vs override)', 'Interpreter','none');
    legend(["follow","override"], 'Interpreter','none', 'Location','best');

    saveas(f, pathPng);
    close(f);
end
