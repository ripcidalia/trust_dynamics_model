function stepA10_behavior_fit_by_participant(run_id, varargin)
% stepA10_behavior_fit_by_participant  Per-participant behavioral model fitting (VALID only).
%
% Step A10 — Per-participant behavioral model fitting
%   - Fits Model 0/1/2 separately for each VALID participant using all doors (≈60).
%   - Computes NLL + Brier, plus AIC/BIC complexity penalties.
%   - Bootstrap (door-resampling with replacement) to quantify parameter uncertainty.
%   - Selects best model per participant via BIC + parsimony rule (ΔBIC <= 2 -> simpler).
%
% This step is mechanistic analysis (not generalization).
%
% Inputs
%   A7 VALID dataset:
%     derived/analysis_runs/<run_id>/stepA7_behavior_dataset/behavior_dataset_valid.mat
%     containing table T with at least:
%       participant_id, door_index, block_index,
%       tau_decision, self_confidence, sc_centered, margin_treshold,
%       followed, is_valid_label
%
% Outputs
%   Writes to:
%     derived/analysis_runs/<run_id>/stepA10_behavior_fit_by_participant/
%       - A10_params_by_participant.csv
%       - A10_params_by_participant.mat   (Tsum + details + meta)
%       - meta.mat / meta.json
%       - figures/*.png
%
% Name-value args
%   "OutDir" (default derived/.../stepA10_behavior_fit_by_participant)
%   "Overwrite" (false)
%   "BootstrapN" (300)
%   "EpsMax" (0.5)
%   "RandomSeed" (1)
%   "DeltaBIC_Parsimony" (2)
%   "KHuge" (100)
%   "MinN" (40)   % safety: if participant has fewer valid doors, still fit but flags it
%
% Dependencies (assumed on path as per your pipeline)
%   must_exist_file, ensure_dir, save_json

    if nargin < 1 || isempty(run_id)
        error("stepA10_behavior_fit_by_participant: run_id is required.");
    end
    run_id = string(run_id);

    p = inputParser;
    p.addParameter("OutDir", "", @(s) isstring(s) || ischar(s));
    p.addParameter("Overwrite", false, @(x) islogical(x) && isscalar(x));
    p.addParameter("BootstrapN", 1000, @(x) isnumeric(x) && isscalar(x) && x>=0);
    p.addParameter("EpsMax", 0.5, @(x) isnumeric(x) && isscalar(x) && x>0 && x<=1);
    p.addParameter("RandomSeed", 1, @(x) isnumeric(x) && isscalar(x));
    p.addParameter("DeltaBIC_Parsimony", 2, @(x) isnumeric(x) && isscalar(x) && x>=0);
    p.addParameter("KHuge", 100, @(x) isnumeric(x) && isscalar(x) && x>0);
    p.addParameter("MinN", 40, @(x) isnumeric(x) && isscalar(x) && x>=1);
    p.parse(varargin{:});
    args = p.Results;

    rng(args.RandomSeed);

    % ------------------------------------------------------------
    % Load A7 VALID dataset
    % ------------------------------------------------------------
    a7Dir = fullfile("derived","analysis_runs",run_id,"stepA7_behavior_dataset");
    validMat = fullfile(a7Dir, "behavior_dataset_valid.mat");
    must_exist_file(validMat, "A7 VALID dataset");

    S = load(validMat, "T");
    if ~isfield(S,"T") || ~istable(S.T)
        error("[A10] VALID mat missing table T.");
    end
    T = S.T;

    reqCols = ["participant_id","door_index","block_index", ...
               "tau_decision","self_confidence","sc_centered","margin_treshold", ...
               "followed","is_valid_label"];
    assert(all(ismember(reqCols, string(T.Properties.VariableNames))), "[A10] A7 VALID missing required columns.");

    % Filter to valid labels
    T = T(T.is_valid_label==1, :);

    % Drop NaNs in predictors used by models
    okPred = isfinite(T.tau_decision) & isfinite(T.self_confidence) & isfinite(T.sc_centered) & isfinite(T.margin_treshold);
    T = T(okPred, :);

    % ------------------------------------------------------------
    % Output directories
    % ------------------------------------------------------------
    outDir = string(args.OutDir);
    if strlength(outDir)==0
        outDir = fullfile("derived","analysis_runs",run_id,"stepA10_behavior_fit_by_participant");
    end
    ensure_dir(outDir);
    figDir = fullfile(outDir, "figures");
    ensure_dir(figDir);

    outCsv  = fullfile(outDir, "A10_params_by_participant.csv");
    outMat  = fullfile(outDir, "A10_params_by_participant.mat");
    metaMat = fullfile(outDir, "meta.mat");
    metaJson= fullfile(outDir, "meta.json");

    if ~args.Overwrite
        if isfile(outCsv) || isfile(outMat) || isfile(metaMat)
            error("[A10] Outputs exist. Set Overwrite=true to replace. (%s)", outDir);
        end
    end

    % ------------------------------------------------------------
    % Per-participant loop (preserve door order in stored sequences)
    % ------------------------------------------------------------
    pid_all = string(T.participant_id);
    uniqP = unique(pid_all);
    nP = numel(uniqP);

    details = struct();
    details.run_id = char(run_id);
    details.participants = cell(nP,1);

    % Summary table (one row per participant)
    Tsum = table();
    Tsum.participant_id = uniqP;
    Tsum.N = zeros(nP,1);

    % Model 0 metrics
    Tsum.m0_NLL = NaN(nP,1);
    Tsum.m0_Brier = NaN(nP,1);
    Tsum.m0_AIC = NaN(nP,1);
    Tsum.m0_BIC = NaN(nP,1);

    % Model 1 params + metrics + CI + flags
    Tsum.m1_k_hat = NaN(nP,1);
    Tsum.m1_NLL = NaN(nP,1);
    Tsum.m1_Brier = NaN(nP,1);
    Tsum.m1_AIC = NaN(nP,1);
    Tsum.m1_BIC = NaN(nP,1);
    Tsum.m1_k_ci_lo = NaN(nP,1);
    Tsum.m1_k_ci_hi = NaN(nP,1);
    Tsum.m1_boot_fail_rate = NaN(nP,1);
    Tsum.m1_flag_k_huge = false(nP,1);
    Tsum.m1_flag_ci_wide = false(nP,1);

    % Model 2 params + metrics + CI + flags
    Tsum.m2_k_hat = NaN(nP,1);
    Tsum.m2_beta_hat = NaN(nP,1);
    Tsum.m2_eps_hat = NaN(nP,1);
    Tsum.m2_NLL = NaN(nP,1);
    Tsum.m2_Brier = NaN(nP,1);
    Tsum.m2_AIC = NaN(nP,1);
    Tsum.m2_BIC = NaN(nP,1);
    Tsum.m2_k_ci_lo = NaN(nP,1);
    Tsum.m2_k_ci_hi = NaN(nP,1);
    Tsum.m2_beta_ci_lo = NaN(nP,1);
    Tsum.m2_beta_ci_hi = NaN(nP,1);
    Tsum.m2_eps_ci_lo = NaN(nP,1);
    Tsum.m2_eps_ci_hi = NaN(nP,1);
    Tsum.m2_boot_fail_rate = NaN(nP,1);
    Tsum.m2_flag_k_huge = false(nP,1);
    Tsum.m2_flag_eps_at_bound = false(nP,1);
    Tsum.m2_flag_ci_wide = false(nP,1);

    % Best model selection outputs
    Tsum.best_model = strings(nP,1);
    Tsum.best_model_idx = NaN(nP,1);
    Tsum.best_model_BIC = NaN(nP,1);
    Tsum.best_model_is_suspect = false(nP,1);
    Tsum.deltaBIC_m1_minus_m0 = NaN(nP,1);
    Tsum.deltaBIC_m2_minus_m1 = NaN(nP,1);
    Tsum.flag_low_N = false(nP,1);

    epsMax = double(args.EpsMax);
    B = double(args.BootstrapN);

    % thresholds for "wide CI" (simple diagnostics)
    CIW_K = 50;       % absolute width
    CIW_BETA = 10;
    CIW_EPS = 0.3;

    for i = 1:nP
        pid = uniqP(i);
        mask = (pid_all == pid);
        Tp = T(mask, :);

        % Preserve temporal order by door_index
        [~,ord] = sort(double(Tp.door_index(:)));
        Tp = Tp(ord,:);

        y = double(Tp.followed(:));          % 1=follow, 0=override
        tau = double(Tp.tau_decision(:));
        m   = double(Tp.margin_treshold(:)); % tau - sc
        scC = double(Tp.sc_centered(:));     % sc - 0.5

        N = numel(y);
        Tsum.N(i) = N;
        if N < args.MinN
            Tsum.flag_low_N(i) = true;
        end

        % ---- Model 0 ----
        p0 = clamp01(tau);
        nll0 = bernoulli_nll(y, p0);
        brier0 = mean((p0 - y).^2);

        Tsum.m0_NLL(i) = nll0;
        Tsum.m0_Brier(i) = brier0;
        Tsum.m0_AIC(i) = aic_from_nll(nll0, 0);
        Tsum.m0_BIC(i) = bic_from_nll(nll0, 0, N);

        % ---- Model 1 (MLE) ----
        % k >= 0 via k = exp(z)
        z0_1 = log(10);
        f1 = @(z) nll_model1(exp(z), m, y);
        [k1_hat, nll1, ok1] = safe_fit_1d_exp(f1, z0_1);

        p1 = sigmoid(k1_hat .* m);
        brier1 = mean((p1 - y).^2);

        Tsum.m1_k_hat(i) = k1_hat;
        Tsum.m1_NLL(i) = nll1;
        Tsum.m1_Brier(i) = brier1;
        Tsum.m1_AIC(i) = aic_from_nll(nll1, 1);
        Tsum.m1_BIC(i) = bic_from_nll(nll1, 1, N);

        % ---- Model 2 (MLE) ----
        % k>=0 via exp(z1), eps in [0,epsMax] via epsMax*sigmoid(z3)
        z0_2 = [log(10); 0; logit(0.05/epsMax)];
        f2 = @(z) nll_model2(exp(z(1)), z(2), epsMax*sigmoid(z(3)), tau, scC, y);
        [k2_hat, beta_hat, eps_hat, nll2, ok2] = safe_fit_model2(f2, z0_2, epsMax);

        z2 = k2_hat.*tau + beta_hat.*scC;
        p2 = (1-eps_hat).*sigmoid(z2) + eps_hat.*0.5;
        p2 = clamp01(p2);
        brier2 = mean((p2 - y).^2);

        Tsum.m2_k_hat(i) = k2_hat;
        Tsum.m2_beta_hat(i) = beta_hat;
        Tsum.m2_eps_hat(i) = eps_hat;
        Tsum.m2_NLL(i) = nll2;
        Tsum.m2_Brier(i) = brier2;
        Tsum.m2_AIC(i) = aic_from_nll(nll2, 3);
        Tsum.m2_BIC(i) = bic_from_nll(nll2, 3, N);

        % ---- Bootstrap uncertainty (door-resampling) ----
        boot = struct();
        boot.B = B;

        if B > 0
            boot.m1_k = NaN(B,1);
            boot.m2_k = NaN(B,1);
            boot.m2_beta = NaN(B,1);
            boot.m2_eps = NaN(B,1);

            fail1 = 0;
            fail2 = 0;

            for b = 1:B
                idx = randi(N, N, 1);  % resample doors with replacement

                yb   = y(idx);
                mb   = m(idx);
                taub = tau(idx);
                scCb = scC(idx);

                % Model 1 bootstrap fit
                f1b = @(z) nll_model1(exp(z), mb, yb);
                [k1b, ~, okb1] = safe_fit_1d_exp(f1b, z0_1);
                if okb1
                    boot.m1_k(b) = k1b;
                else
                    fail1 = fail1 + 1;
                end

                % Model 2 bootstrap fit
                f2b = @(z) nll_model2(exp(z(1)), z(2), epsMax*sigmoid(z(3)), taub, scCb, yb);
                [k2b, betab, epsb, ~, okb2] = safe_fit_model2(f2b, z0_2, epsMax);
                if okb2
                    boot.m2_k(b) = k2b;
                    boot.m2_beta(b) = betab;
                    boot.m2_eps(b) = epsb;
                else
                    fail2 = fail2 + 1;
                end
            end

            boot.m1_fail_rate = fail1 / max(1,B);
            boot.m2_fail_rate = fail2 / max(1,B);

            Tsum.m1_boot_fail_rate(i) = boot.m1_fail_rate;
            Tsum.m2_boot_fail_rate(i) = boot.m2_fail_rate;

            % CI (percentile) ignoring NaNs
            [lo,hi] = ci95_percentile(boot.m1_k);
            Tsum.m1_k_ci_lo(i) = lo;
            Tsum.m1_k_ci_hi(i) = hi;

            [lo,hi] = ci95_percentile(boot.m2_k);
            Tsum.m2_k_ci_lo(i) = lo;
            Tsum.m2_k_ci_hi(i) = hi;

            [lo,hi] = ci95_percentile(boot.m2_beta);
            Tsum.m2_beta_ci_lo(i) = lo;
            Tsum.m2_beta_ci_hi(i) = hi;

            [lo,hi] = ci95_percentile(boot.m2_eps);
            Tsum.m2_eps_ci_lo(i) = lo;
            Tsum.m2_eps_ci_hi(i) = hi;
        end

        % ---- Flags (simple diagnostics) ----
        Tsum.m1_flag_k_huge(i) = isfinite(k1_hat) && (k1_hat > args.KHuge);
        Tsum.m2_flag_k_huge(i) = isfinite(k2_hat) && (k2_hat > args.KHuge);
        Tsum.m2_flag_eps_at_bound(i) = isfinite(eps_hat) && (eps_hat < 1e-3 || eps_hat > (epsMax - 1e-3));

        % CI width flags (only if CI finite)
        if isfinite(Tsum.m1_k_ci_lo(i)) && isfinite(Tsum.m1_k_ci_hi(i))
            Tsum.m1_flag_ci_wide(i) = (Tsum.m1_k_ci_hi(i) - Tsum.m1_k_ci_lo(i)) > CIW_K;
        else
            Tsum.m1_flag_ci_wide(i) = true;
        end

        wide2 = false;
        if isfinite(Tsum.m2_k_ci_lo(i)) && isfinite(Tsum.m2_k_ci_hi(i))
            wide2 = wide2 || ((Tsum.m2_k_ci_hi(i) - Tsum.m2_k_ci_lo(i)) > CIW_K);
        else
            wide2 = true;
        end
        if isfinite(Tsum.m2_beta_ci_lo(i)) && isfinite(Tsum.m2_beta_ci_hi(i))
            wide2 = wide2 || ((Tsum.m2_beta_ci_hi(i) - Tsum.m2_beta_ci_lo(i)) > CIW_BETA);
        else
            wide2 = true;
        end
        if isfinite(Tsum.m2_eps_ci_lo(i)) && isfinite(Tsum.m2_eps_ci_hi(i))
            wide2 = wide2 || ((Tsum.m2_eps_ci_hi(i) - Tsum.m2_eps_ci_lo(i)) > CIW_EPS);
        else
            wide2 = true;
        end
        Tsum.m2_flag_ci_wide(i) = wide2;

        % Additional suspect signals: optimizer failure or high bootstrap fail rate
        suspect1 = (~ok1) || (isfinite(Tsum.m1_boot_fail_rate(i)) && Tsum.m1_boot_fail_rate(i) > 0.2) || Tsum.m1_flag_ci_wide(i) || Tsum.m1_flag_k_huge(i);
        suspect2 = (~ok2) || (isfinite(Tsum.m2_boot_fail_rate(i)) && Tsum.m2_boot_fail_rate(i) > 0.2) || Tsum.m2_flag_ci_wide(i) || Tsum.m2_flag_k_huge(i) || Tsum.m2_flag_eps_at_bound(i);

        % ---- Best model selection (BIC + parsimony) ----
        bic0 = Tsum.m0_BIC(i);
        bic1 = Tsum.m1_BIC(i);
        bic2 = Tsum.m2_BIC(i);

        Tsum.deltaBIC_m1_minus_m0(i) = bic1 - bic0;
        Tsum.deltaBIC_m2_minus_m1(i) = bic2 - bic1;

        [bestIdx, bestName] = select_best_model_bic_parsimony([bic0 bic1 bic2], args.DeltaBIC_Parsimony);
        Tsum.best_model_idx(i) = bestIdx;
        Tsum.best_model(i) = bestName;
        bicVec = [bic0 bic1 bic2];
        Tsum.best_model_BIC(i) = bicVec(bestIdx);

        % Best-model suspect flag
        switch bestIdx
            case 1, Tsum.best_model_is_suspect(i) = false;          % Model 0 has no params
            case 2, Tsum.best_model_is_suspect(i) = suspect1;
            case 3, Tsum.best_model_is_suspect(i) = suspect2;
        end

        % ---- Store per-participant details (including sequences) ----
        D = struct();
        D.participant_id = char(pid);
        D.N = N;
        D.seq = struct( ...
            "door_index", double(Tp.door_index(:)), ...
            "block_index", double(Tp.block_index(:)), ...
            "y_follow", y(:), ...
            "tau", tau(:), ...
            "margin_treshold", m(:), ...
            "sc_centered", scC(:));

        D.model0 = struct("nll",nll0,"brier",brier0);
        D.model1 = struct("k_hat",k1_hat,"nll",nll1,"brier",brier1, ...
                          "k_ci95",[Tsum.m1_k_ci_lo(i) Tsum.m1_k_ci_hi(i)], ...
                          "boot_fail_rate",Tsum.m1_boot_fail_rate(i));
        D.model2 = struct("k_hat",k2_hat,"beta_hat",beta_hat,"eps_hat",eps_hat,"nll",nll2,"brier",brier2, ...
                          "k_ci95",[Tsum.m2_k_ci_lo(i) Tsum.m2_k_ci_hi(i)], ...
                          "beta_ci95",[Tsum.m2_beta_ci_lo(i) Tsum.m2_beta_ci_hi(i)], ...
                          "eps_ci95",[Tsum.m2_eps_ci_lo(i) Tsum.m2_eps_ci_hi(i)], ...
                          "boot_fail_rate",Tsum.m2_boot_fail_rate(i));

        D.flags = struct( ...
            "low_N", Tsum.flag_low_N(i), ...
            "m1_k_huge", Tsum.m1_flag_k_huge(i), ...
            "m1_ci_wide", Tsum.m1_flag_ci_wide(i), ...
            "m2_k_huge", Tsum.m2_flag_k_huge(i), ...
            "m2_eps_at_bound", Tsum.m2_flag_eps_at_bound(i), ...
            "m2_ci_wide", Tsum.m2_flag_ci_wide(i), ...
            "best_model_is_suspect", Tsum.best_model_is_suspect(i));

        D.best = struct("model_idx",bestIdx,"model_name",char(bestName),"bic",Tsum.best_model_BIC(i));

        if B > 0
            D.bootstrap = boot; % includes arrays (k/beta/eps) -> can be large but manageable (nP~10)
        end

        details.participants{i} = D;
    end

    % ------------------------------------------------------------
    % Save artifacts
    % ------------------------------------------------------------
    writetable(Tsum, outCsv);
    save(outMat, "Tsum", "details", "-v7.3");

    meta = struct();
    meta.run_id = char(run_id);
    meta.created = char(datetime('now','Format','yyyy-MM-dd HH:mm:ss'));
    meta.valid_participants = nP;
    meta.bootstrapN = B;
    meta.epsMax = epsMax;
    meta.random_seed = args.RandomSeed;
    meta.deltaBIC_parsimony = args.DeltaBIC_Parsimony;
    meta.k_huge_threshold = args.KHuge;
    meta.minN_flag_threshold = args.MinN;
    save(metaMat, "meta");
    save_json(metaJson, meta);

    % ------------------------------------------------------------
    % Figures (minimal, thesis-friendly)
    % ------------------------------------------------------------
    make_fig_bic_by_participant(fullfile(figDir, "bic_by_participant.png"), Tsum);
    make_fig_best_model_counts(fullfile(figDir, "best_model_counts.png"), Tsum);

    fprintf("[Step A10] Done.\n");
    fprintf("  VALID participants: %d\n", nP);
    fprintf("  BootstrapN: %d\n", B);
    fprintf("  Output dir: %s\n", outDir);
    fprintf("  Wrote: %s\n", outCsv);
end

% ======================================================================
% Selection rule: BIC + parsimony (ΔBIC <= delta => simpler)
% Models are ordered by simplicity: 0 < 1 < 2
% ======================================================================
function [bestIdx, bestName] = select_best_model_bic_parsimony(bicVec012, delta)
    % bicVec012: [bic0 bic1 bic2]
    bic0 = bicVec012(1);
    bic1 = bicVec012(2);
    bic2 = bicVec012(3);

    % Start: best by raw BIC
    [~, idxMin] = min([bic0 bic1 bic2]);
    bestIdx = idxMin;

    % Apply parsimony: if within delta of a simpler model, choose simpler
    % Check simplicity chain 0 -> 1 -> 2
    % If model2 wins but within delta of model1, choose model1.
    if bestIdx == 3 && isfinite(bic2) && isfinite(bic1) && (bic2 - bic1) <= delta
        bestIdx = 2;
    end
    % If model1 wins but within delta of model0, choose model0.
    if bestIdx == 2 && isfinite(bic1) && isfinite(bic0) && (bic1 - bic0) <= delta
        bestIdx = 1;
    end
    % If model2 wins but within delta of model0 directly (rare), still choose 0.
    if bestIdx == 3 && isfinite(bic2) && isfinite(bic0) && (bic2 - bic0) <= delta
        bestIdx = 1;
    end

    switch bestIdx
        case 1, bestName = "model0_trust_as_probability";
        case 2, bestName = "model1_threshold";
        case 3, bestName = "model2_offset_lapse";
        otherwise, bestName = "unknown";
    end
end

% ======================================================================
% Model likelihoods (same structure as A8; adapted for per-participant)
% ======================================================================
function nll = nll_model1(k, margin, y_follow)
    if ~isfinite(k) || k < 0
        nll = Inf; return;
    end
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

% ======================================================================
% Fitting wrappers (safe guards)
% ======================================================================
function [k_hat, nll_hat, ok] = safe_fit_1d_exp(f, z0)
    ok = true;
    try
        zhat = fminsearch(f, z0, optimset('Display','off'));
        k_hat = exp(zhat);
        nll_hat = f(zhat);
        if ~isfinite(k_hat) || ~isfinite(nll_hat)
            ok = false;
        end
    catch
        ok = false;
        k_hat = NaN;
        nll_hat = Inf;
    end
end

function [k_hat, beta_hat, eps_hat, nll_hat, ok] = safe_fit_model2(f, z0, epsMax)
    ok = true;
    try
        zhat = fminsearch(f, z0, optimset('Display','off'));
        k_hat    = exp(zhat(1));
        beta_hat = zhat(2);
        eps_hat  = epsMax * sigmoid(zhat(3));
        nll_hat  = f(zhat);
        if ~isfinite(k_hat) || ~isfinite(beta_hat) || ~isfinite(eps_hat) || ~isfinite(nll_hat)
            ok = false;
        end
    catch
        ok = false;
        k_hat = NaN; beta_hat = NaN; eps_hat = NaN; nll_hat = Inf;
    end
end

% ======================================================================
% Metrics helpers
% ======================================================================
function a = aic_from_nll(nll, k_params)
    a = 2*k_params + 2*nll;
end

function b = bic_from_nll(nll, k_params, N)
    b = k_params*log(max(1,N)) + 2*nll;
end

function [lo,hi] = ci95_percentile(x)
    x = x(:);
    x = x(isfinite(x));
    if isempty(x)
        lo = NaN; hi = NaN; return;
    end
    qs = prctile(x, [2.5 97.5]);
    lo = qs(1); hi = qs(2);
end

% ======================================================================
% Simple plots (minimal, consistent with A8/A9 style)
% ======================================================================
function make_fig_bic_by_participant(pathPng, Tsum)
    f = figure('Visible','off');

    pid = string(Tsum.participant_id);
    x = 1:numel(pid);

    bic0 = Tsum.m0_BIC;
    bic1 = Tsum.m1_BIC;
    bic2 = Tsum.m2_BIC;

    % Grouped bars
    M = [bic0 bic1 bic2];
    bar(M);
    grid on;
    xlabel('Participant');
    ylabel('BIC (lower is better)');
    title('A10: BIC by participant (Models 0/1/2)', 'Interpreter','none');

    % Use participant labels (rotated)
    set(gca,'XTick',x,'XTickLabel',pid);
    xtickangle(30);

    legend({'Model0','Model1','Model2'}, 'Location','best');
    saveas(f, pathPng);
    close(f);
end

function make_fig_best_model_counts(pathPng, Tsum)
    f = figure('Visible','off');
    grid on; hold on;

    names = ["model0_trust_as_probability","model1_threshold","model2_offset_lapse"];
    counts = zeros(1,numel(names));
    for i = 1:numel(names)
        counts(i) = sum(string(Tsum.best_model)==names(i));
    end

    bar(counts);
    set(gca,'XTick',1:numel(names),'XTickLabel',{'Model0','Model1','Model2'});
    ylabel('Count');
    title('A10: Best model counts across VALID participants (BIC + parsimony)', 'Interpreter','none');
    saveas(f, pathPng);
    close(f);
end

% ======================================================================
% Math helpers
% ======================================================================
function p = clamp01(x)
    p = min(max(x, 0), 1);
end

function s = sigmoid(z)
    s = 1 ./ (1 + exp(-z));
end

function z = logit(p)
    p = min(max(p, 1e-12), 1-1e-12);
    z = log(p./(1-p));
end
