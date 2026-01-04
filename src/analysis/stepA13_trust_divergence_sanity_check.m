function stepA13_trust_divergence_sanity_check(run_id, varargin)
% stepA13_trust_divergence_sanity_check  Trust realism sanity check + fallback audit.
%
% Step A13 — Trust realism sanity check (continuous-time divergence on full time grid)
%
% Compares trust trajectories on the FULL simulation grid (not only door times):
%
%   A) simple-mode replay trust vs coupled-mode trust (global params and personalized)
%   B) coupled trust under GLOBAL behavior params vs coupled trust under PERSONALIZED params
%      (with MATCHED random seeds to isolate parameter effects from sampling noise)
%
% Key divergence metrics (per rollout):
%   e(t_k) = tau_cpl(t_k) - tau_simple(t_k)   [option A]
%   e(t_k) = tau_cpl_personal(t_k) - tau_cpl_global(t_k) [option B]
%
%   IAD   = sum_k |e(t_k)| * dt                          [trust*s]
%   MAE_t = (1/T) * sum_k |e(t_k)| * dt                  [trust]
%   RMSE_t= sqrt( (1/T) * sum_k e(t_k)^2 * dt )          [trust]
%   MaxAbs= max_k |e(t_k)|                               [trust]
%
% Outputs (writes)
%   derived/analysis_runs/<run_id>/stepA13_trust_divergence_sanity/
%     - A13_divergence_by_participant.csv
%     - A13_divergence_by_participant.mat   (T + meta + pooled curve bundles)
%     - A13_rollout_level_metrics.mat       (optional, large but useful)
%     - meta.mat / meta.json
%     - figures/*.png
%     - figures/by_participant/*.png
%
% Inputs (reads)
%   - A1 / A3 for theta_star, dt, VALID participants (loaded like A9/A12)
%   - A8 global behavior fit:
%       derived/.../stepA8_behavior_fit_eval/fit_params.mat
%   - A10 participant best models/params:
%       derived/.../stepA10_behavior_fit_by_participant/A10_params_by_participant.mat
%   - A11 robustness (optional guardrails):
%       derived/.../stepA11_behavior_param_robustness/A11_blockwise_params.mat
%
% Name-value args
%   "OutDir" (default derived/.../stepA13_trust_divergence_sanity)
%   "Overwrite" (false)
%   "RolloutsPerParticipant" (1000)
%   "Quantiles" ([0.05 0.95])
%   "RandomSeed" (1)
%   "UseA11Guards" (true)
%   "FallbackStrategy" ("simple")           % 2->1->0
%   "EpsMax" (0.5)
%   "GlobalModelIdx" (2)                    % compare against global Model 2 by default
%   "PooledTimeGridN" (201)                 % normalized-time points for pooled curves
%   "SaveRolloutLevelMat" (true)            % can be big
%
% Notes
%   - This step re-simulates trajectories to capture full tau_hist(t).
%   - For pooled plots across participants with different time horizons,
%     we resample error curves onto a normalized-time grid s in [0,1].
%
% Dependencies (assumed on path)
%   must_exist_file, ensure_dir, save_json, load_participants_struct, find_theta_in_struct
%   trust_simulate_or_predict_one_participant, behavioral_model
%
    if nargin < 1 || isempty(run_id)
        error("stepA13_trust_divergence_sanity_check: run_id is required.");
    end
    run_id = string(run_id);

    p = inputParser;
    p.addParameter("OutDir", "", @(s) isstring(s) || ischar(s));
    p.addParameter("Overwrite", false, @(x) islogical(x) && isscalar(x));
    p.addParameter("RolloutsPerParticipant", 1000, @(x) isnumeric(x) && isscalar(x) && x>=1);
    p.addParameter("Quantiles", [0.05 0.95], @(x) isnumeric(x) && numel(x)==2 && all(x>=0) && all(x<=1));
    p.addParameter("RandomSeed", 1, @(x) isnumeric(x) && isscalar(x));
    p.addParameter("UseA11Guards", true, @(x) islogical(x) && isscalar(x));
    p.addParameter("FallbackStrategy", "simple", @(s) isstring(s) || ischar(s));
    p.addParameter("EpsMax", 0.5, @(x) isnumeric(x) && isscalar(x) && x>0 && x<=1);
    p.addParameter("GlobalModelIdx", 2, @(x) isnumeric(x) && isscalar(x) && ismember(round(x), [0 1 2]));
    p.addParameter("PooledTimeGridN", 201, @(x) isnumeric(x) && isscalar(x) && x>=51);
    p.addParameter("SaveRolloutLevelMat", true, @(x) islogical(x) && isscalar(x));
    p.parse(varargin{:});
    args = p.Results;

    rng(args.RandomSeed);

    % ------------------------------------------------------------
    % Locate inputs
    % ------------------------------------------------------------
    a8Dir = fullfile("derived","analysis_runs",run_id,"stepA8_behavior_fit_eval");
    fitMatA8 = fullfile(a8Dir, "fit_params.mat");
    must_exist_file(fitMatA8, "A8 fit_params.mat (global behavior params)");

    a10Dir = fullfile("derived","analysis_runs",run_id,"stepA10_behavior_fit_by_participant");
    a10Mat = fullfile(a10Dir, "A10_params_by_participant.mat");
    must_exist_file(a10Mat, "A10_params_by_participant.mat");

    a11Dir = fullfile("derived","analysis_runs",run_id,"stepA11_behavior_param_robustness");
    a11Mat = fullfile(a11Dir, "A11_blockwise_params.mat");
    haveA11 = isfile(a11Mat);

    % ------------------------------------------------------------
    % Load theta_star, dt, VALID participants (like A9/A12)
    % ------------------------------------------------------------
    [theta_star, dt, validParticipants] = local_load_theta_dt_and_valid_participants_like_A5(run_id);

    % ------------------------------------------------------------
    % Load A8 global fit
    % ------------------------------------------------------------
    Sfit = load(fitMatA8, "fit");
    if ~isfield(Sfit,"fit") || ~isstruct(Sfit.fit)
        error("[A13] A8 fit_params.mat missing 'fit' struct.");
    end
    fit = Sfit.fit;

    globalModelIdx = round(args.GlobalModelIdx);
    [globalName, bpar_global] = resolve_behavior_params_global(globalModelIdx, fit);

    % ------------------------------------------------------------
    % Load A10 table (best model + params)
    % ------------------------------------------------------------
    Sa10 = load(a10Mat);
    Ta10 = local_find_first_table(Sa10, ["Tsum","T"]);
    if isempty(Ta10)
        error("[A13] Could not find A10 summary table in %s (expected Tsum).", a10Mat);
    end
    if ~ismember("participant_id", string(Ta10.Properties.VariableNames)) || ~ismember("best_model_idx", string(Ta10.Properties.VariableNames))
        error("[A13] A10 table must contain participant_id and best_model_idx.");
    end

    % ------------------------------------------------------------
    % Load A11 (optional)
    % ------------------------------------------------------------
    Tblocks = table();
    if haveA11
        Sa11 = load(a11Mat);
        Tblocks = local_find_first_table(Sa11, ["Tblocks","T"]);
        if isempty(Tblocks)
            Tblocks = table(); % ok
        end
    end

    % ------------------------------------------------------------
    % Output directories
    % ------------------------------------------------------------
    outDir = string(args.OutDir);
    if strlength(outDir)==0
        outDir = fullfile("derived","analysis_runs",run_id,"stepA13_trust_divergence_sanity");
    end
    ensure_dir(outDir);
    figDir = fullfile(outDir, "figures"); ensure_dir(figDir);
    figByP = fullfile(figDir, "by_participant"); ensure_dir(figByP);

    outCsv = fullfile(outDir, "A13_divergence_by_participant.csv");
    outMat = fullfile(outDir, "A13_divergence_by_participant.mat");
    rolloutMat = fullfile(outDir, "A13_rollout_level_metrics.mat");
    metaMat = fullfile(outDir, "meta.mat");
    metaJson= fullfile(outDir, "meta.json");

    if ~args.Overwrite
        if isfile(outCsv) || isfile(outMat) || (args.SaveRolloutLevelMat && isfile(rolloutMat))
            error("[A13] Outputs exist. Set Overwrite=true to replace. (%s)", outDir);
        end
    end

    % ------------------------------------------------------------
    % Participants
    % ------------------------------------------------------------
    uniqP = local_get_participant_ids(validParticipants);
    nP = numel(uniqP);

    R   = double(args.RolloutsPerParticipant);
    qlo = double(args.Quantiles(1));
    qhi = double(args.Quantiles(2));

    fprintf("[A13] Trust divergence sanity check\n");
    fprintf("      VALID participants: %d\n", nP);
    fprintf("      Rollouts/participant: %d\n", R);
    fprintf("      dt: %.6g s\n", dt);
    fprintf("      Global behavior: idx=%d (%s)\n", globalModelIdx, globalName);

    % ------------------------------------------------------------
    % Allocate participant summary table
    % ------------------------------------------------------------
    T = table();
    T.participant_id = uniqP;

    % A) Simple vs Coupled (GLOBAL) metrics summary
    T.A_simple_vs_cpl_global_MAE_mean  = NaN(nP,1);
    T.A_simple_vs_cpl_global_MAE_qlo   = NaN(nP,1);
    T.A_simple_vs_cpl_global_MAE_qhi   = NaN(nP,1);
    T.A_simple_vs_cpl_global_RMSE_mean = NaN(nP,1);
    T.A_simple_vs_cpl_global_RMSE_qlo  = NaN(nP,1);
    T.A_simple_vs_cpl_global_RMSE_qhi  = NaN(nP,1);
    T.A_simple_vs_cpl_global_IAD_mean  = NaN(nP,1);
    T.A_simple_vs_cpl_global_IAD_qlo   = NaN(nP,1);
    T.A_simple_vs_cpl_global_IAD_qhi   = NaN(nP,1);
    T.A_simple_vs_cpl_global_MaxAbs_mean = NaN(nP,1);
    T.A_simple_vs_cpl_global_MaxAbs_qlo  = NaN(nP,1);
    T.A_simple_vs_cpl_global_MaxAbs_qhi  = NaN(nP,1);

    % A) Simple vs Coupled (PERSONALIZED) metrics summary
    T.A_simple_vs_cpl_person_MAE_mean  = NaN(nP,1);
    T.A_simple_vs_cpl_person_MAE_qlo   = NaN(nP,1);
    T.A_simple_vs_cpl_person_MAE_qhi   = NaN(nP,1);
    T.A_simple_vs_cpl_person_RMSE_mean = NaN(nP,1);
    T.A_simple_vs_cpl_person_RMSE_qlo  = NaN(nP,1);
    T.A_simple_vs_cpl_person_RMSE_qhi  = NaN(nP,1);
    T.A_simple_vs_cpl_person_IAD_mean  = NaN(nP,1);
    T.A_simple_vs_cpl_person_IAD_qlo   = NaN(nP,1);
    T.A_simple_vs_cpl_person_IAD_qhi   = NaN(nP,1);
    T.A_simple_vs_cpl_person_MaxAbs_mean = NaN(nP,1);
    T.A_simple_vs_cpl_person_MaxAbs_qlo  = NaN(nP,1);
    T.A_simple_vs_cpl_person_MaxAbs_qhi  = NaN(nP,1);

    % B) Coupled PERSONALIZED vs Coupled GLOBAL (matched seeds)
    T.B_cpl_person_vs_global_MAE_mean  = NaN(nP,1);
    T.B_cpl_person_vs_global_MAE_qlo   = NaN(nP,1);
    T.B_cpl_person_vs_global_MAE_qhi   = NaN(nP,1);
    T.B_cpl_person_vs_global_RMSE_mean = NaN(nP,1);
    T.B_cpl_person_vs_global_RMSE_qlo  = NaN(nP,1);
    T.B_cpl_person_vs_global_RMSE_qhi  = NaN(nP,1);
    T.B_cpl_person_vs_global_IAD_mean  = NaN(nP,1);
    T.B_cpl_person_vs_global_IAD_qlo   = NaN(nP,1);
    T.B_cpl_person_vs_global_IAD_qhi   = NaN(nP,1);
    T.B_cpl_person_vs_global_MaxAbs_mean = NaN(nP,1);
    T.B_cpl_person_vs_global_MaxAbs_qlo  = NaN(nP,1);
    T.B_cpl_person_vs_global_MaxAbs_qhi  = NaN(nP,1);

    % Guardrail/fallback audit (personalized resolver)
    T.model_idx_a10 = NaN(nP,1);
    T.model_idx_used = NaN(nP,1);
    T.fallback_applied = false(nP,1);
    T.fallback_reason = strings(nP,1);
    T.model_name_used = strings(nP,1);

    % ------------------------------------------------------------
    % Pooled curve storage on normalized time grid s in [0,1]
    %   For A (simple vs coupled global/personal) and B (personal vs global)
    %   We store participant-level rollout-mean curves, then pool across participants.
    % ------------------------------------------------------------
    Ns = round(args.PooledTimeGridN);
    sgrid = linspace(0,1,Ns)';

    E_A_global_byP  = NaN(Ns, nP);
    E_A_person_byP  = NaN(Ns, nP);
    E_B_byP         = NaN(Ns, nP);

    Acc_A_global_byP = NaN(Ns, nP);
    Acc_A_person_byP = NaN(Ns, nP);
    Acc_B_byP        = NaN(Ns, nP);

    % Optional rollout-level metric arrays for debugging / deeper analysis
    rollLevel = struct();
    rollLevel.meta = struct("run_id",char(run_id),"R",R,"dt",dt,"global_model_idx",globalModelIdx, ...
        "global_model_name",char(globalName),"random_seed",args.RandomSeed);
    if args.SaveRolloutLevelMat
        rollLevel.pid = uniqP;
        rollLevel.A_global = init_rollout_metric_store(nP, R);
        rollLevel.A_person = init_rollout_metric_store(nP, R);
        rollLevel.B = init_rollout_metric_store(nP, R);
    end

    % ------------------------------------------------------------
    % Main loop
    % ------------------------------------------------------------
    for pi = 1:nP
        pid = uniqP(pi);
        Pp = get_participant_from_collection(validParticipants, pid);

        % Personalized params (with guardrail fallback)
        [bpar_personal, info] = resolve_personalized_behavior_params( ...
            pid, Ta10, Tblocks, args.UseA11Guards, args.EpsMax, args.FallbackStrategy);

        T.model_idx_a10(pi) = info.model_idx_a10;
        T.model_idx_used(pi)= info.model_idx_used;
        T.model_name_used(pi)= info.model_name_used;
        T.fallback_applied(pi)= info.fallback_applied;
        T.fallback_reason(pi)= info.fallback_reason;

        % Simple replay (deterministic)
        simSimple = trust_simulate_or_predict_one_participant("simple", theta_star, Pp, dt);
        tauS = double(simSimple.tau_hist(:));
        tgrid = double(simSimple.t_grid(:));
        if isempty(tgrid) || isempty(tauS) || numel(tgrid) ~= numel(tauS)
            error("[A13] pid=%s returned invalid simple trajectory.", pid);
        end
        Ttotal = max(tgrid) - min(tgrid);
        if ~isfinite(Ttotal) || Ttotal <= 0
            % fallback, but should not happen
            Ttotal = (numel(tgrid)-1) * dt;
        end

        % Per-rollout metrics arrays
        Aglob = NaN(R,4); % [MAE RMSE IAD MaxAbs]
        Apers = NaN(R,4);
        Bpg   = NaN(R,4);

        % For pooled curves: use rollout-mean error curve per participant
        eAglob_mean = zeros(numel(tgrid),1);
        eApers_mean = zeros(numel(tgrid),1);
        eB_mean     = zeros(numel(tgrid),1);

        accAglob_mean = zeros(numel(tgrid),1);
        accApers_mean = zeros(numel(tgrid),1);
        accB_mean     = zeros(numel(tgrid),1);

        for r = 1:R
            % Matched seed across global vs personalized coupled sims for B.
            seed = args.RandomSeed + 100000*pi + r;
            % --- Coupled GLOBAL ---
            rng(seed);
            simG = trust_simulate_or_predict_one_participant("coupled", theta_star, Pp, dt, bpar_global);
            tauG = double(simG.tau_hist(:));

            % --- Coupled PERSONALIZED ---
            rng(seed);
            simP = trust_simulate_or_predict_one_participant("coupled", theta_star, Pp, dt, bpar_personal);
            tauP = double(simP.tau_hist(:));

            % Sanity: lengths match simple
            tauG = local_align_length(tauG, tauS);
            tauP = local_align_length(tauP, tauS);

            % A-global: tauG - tauS
            eA_g = tauG - tauS;
            [mae, rmse, iad, mx, acc] = divergence_metrics_from_error(eA_g, dt, Ttotal);
            Aglob(r,:) = [mae rmse iad mx];

            % A-personal: tauP - tauS
            eA_p = tauP - tauS;
            [mae, rmse, iad, mx, accp] = divergence_metrics_from_error(eA_p, dt, Ttotal);
            Apers(r,:) = [mae rmse iad mx];

            % B: tauP - tauG
            eB = tauP - tauG;
            [mae, rmse, iad, mx, accb] = divergence_metrics_from_error(eB, dt, Ttotal);
            Bpg(r,:) = [mae rmse iad mx];

            % rollout-mean curves
            eAglob_mean = eAglob_mean + eA_g;
            eApers_mean = eApers_mean + eA_p;
            eB_mean     = eB_mean     + eB;

            accAglob_mean = accAglob_mean + acc;
            accApers_mean = accApers_mean + accp;
            accB_mean     = accB_mean     + accb;

            % Store rollout-level metrics if requested
            if args.SaveRolloutLevelMat
                rollLevel.A_global.MAE(pi,r)  = Aglob(r,1);
                rollLevel.A_global.RMSE(pi,r) = Aglob(r,2);
                rollLevel.A_global.IAD(pi,r)  = Aglob(r,3);
                rollLevel.A_global.MaxAbs(pi,r)=Aglob(r,4);

                rollLevel.A_person.MAE(pi,r)  = Apers(r,1);
                rollLevel.A_person.RMSE(pi,r) = Apers(r,2);
                rollLevel.A_person.IAD(pi,r)  = Apers(r,3);
                rollLevel.A_person.MaxAbs(pi,r)=Apers(r,4);

                rollLevel.B.MAE(pi,r)  = Bpg(r,1);
                rollLevel.B.RMSE(pi,r) = Bpg(r,2);
                rollLevel.B.IAD(pi,r)  = Bpg(r,3);
                rollLevel.B.MaxAbs(pi,r)=Bpg(r,4);
            end
        end

        % Finish mean curves for this participant
        eAglob_mean = eAglob_mean / R;
        eApers_mean = eApers_mean / R;
        eB_mean     = eB_mean     / R;

        accAglob_mean = accAglob_mean / R;
        accApers_mean = accApers_mean / R;
        accB_mean     = accB_mean     / R;

        % Resample participant mean curves to normalized time sgrid for pooling
        s = (tgrid - min(tgrid)) ./ max(1e-12, (max(tgrid)-min(tgrid)));
        E_A_global_byP(:,pi)  = interp1(s, eAglob_mean, sgrid, "linear", "extrap");
        E_A_person_byP(:,pi)  = interp1(s, eApers_mean, sgrid, "linear", "extrap");
        E_B_byP(:,pi)         = interp1(s, eB_mean,    sgrid, "linear", "extrap");

        % For accumulated curves, resample the accumulated |e| integral (already in trust*s)
        Acc_A_global_byP(:,pi)= interp1(s, accAglob_mean, sgrid, "linear", "extrap");
        Acc_A_person_byP(:,pi)= interp1(s, accApers_mean, sgrid, "linear", "extrap");
        Acc_B_byP(:,pi)       = interp1(s, accB_mean,     sgrid, "linear", "extrap");

        % Participant summaries (quantiles across rollouts)
        T = local_write_summary_row(T, pi, Aglob, "A_simple_vs_cpl_global", qlo, qhi);
        T = local_write_summary_row(T, pi, Apers, "A_simple_vs_cpl_person", qlo, qhi);
        T = local_write_summary_row(T, pi, Bpg,   "B_cpl_person_vs_global", qlo, qhi);

        % Per-participant figures (real time grid)
        make_fig_error_and_accum(fullfile(figByP, sprintf("pid_%s_A_global_error.png", sanitize_pid(pid))), ...
            pid, tgrid, eAglob_mean, accAglob_mean, "A: coupled(GLOBAL) - simple");
        make_fig_error_and_accum(fullfile(figByP, sprintf("pid_%s_A_person_error.png", sanitize_pid(pid))), ...
            pid, tgrid, eApers_mean, accApers_mean, "A: coupled(PERSONAL) - simple");
        make_fig_error_and_accum(fullfile(figByP, sprintf("pid_%s_B_person_minus_global.png", sanitize_pid(pid))), ...
            pid, tgrid, eB_mean, accB_mean, "B: coupled(PERSONAL) - coupled(GLOBAL)");
    end

    % ------------------------------------------------------------
    % Pooled curves + pooled figures (normalized time)
    % ------------------------------------------------------------
    pooled = struct();
    pooled.sgrid = sgrid;

    pooled.A_global = summarize_curves(E_A_global_byP, Acc_A_global_byP, qlo, qhi);
    pooled.A_person = summarize_curves(E_A_person_byP, Acc_A_person_byP, qlo, qhi);
    pooled.B        = summarize_curves(E_B_byP,        Acc_B_byP,        qlo, qhi);

    make_fig_pooled_error(fullfile(figDir, "pooled_A_error_simple_vs_coupled_global.png"), ...
        sgrid, pooled.A_global, "Pooled error e(t): coupled(GLOBAL) - simple (normalized time)");
    make_fig_pooled_error(fullfile(figDir, "pooled_A_error_simple_vs_coupled_personal.png"), ...
        sgrid, pooled.A_person, "Pooled error e(t): coupled(PERSONAL) - simple (normalized time)");
    make_fig_pooled_error(fullfile(figDir, "pooled_B_error_person_minus_global.png"), ...
        sgrid, pooled.B, "Pooled error e(t): coupled(PERSONAL) - coupled(GLOBAL) (normalized time)");

    make_fig_pooled_accum(fullfile(figDir, "pooled_A_accum_abs_error_global.png"), ...
        sgrid, pooled.A_global, "Pooled accumulated |e| dt: coupled(GLOBAL) - simple (trust*s)");
    make_fig_pooled_accum(fullfile(figDir, "pooled_A_accum_abs_error_personal.png"), ...
        sgrid, pooled.A_person, "Pooled accumulated |e| dt: coupled(PERSONAL) - simple (trust*s)");
    make_fig_pooled_accum(fullfile(figDir, "pooled_B_accum_abs_error_person_minus_global.png"), ...
        sgrid, pooled.B, "Pooled accumulated |e| dt: coupled(PERSONAL) - coupled(GLOBAL) (trust*s)");

    % ------------------------------------------------------------
    % Save outputs
    % ------------------------------------------------------------
    writetable(T, outCsv);

    meta = struct();
    meta.run_id = char(run_id);
    meta.created = char(datetime('now','Format','yyyy-MM-dd HH:mm:ss'));
    meta.valid_participants = nP;
    meta.rollouts_per_participant = R;
    meta.quantiles = [qlo qhi];
    meta.random_seed = args.RandomSeed;
    meta.dt = dt;
    meta.theta_dim = numel(theta_star);
    meta.global_model_idx = globalModelIdx;
    meta.global_model_name = char(globalName);
    meta.use_a11_guards = args.UseA11Guards && haveA11;
    meta.have_a11_file = haveA11;
    meta.a8_fit_file = char(fitMatA8);
    meta.a10_file = char(a10Mat);
    meta.a11_file = char(a11Mat);
    meta.fallback_strategy = char(string(args.FallbackStrategy));
    meta.pooled_time_grid_n = Ns;

    save(outMat, "T", "meta", "pooled", "-v7.3");
    save(metaMat, "meta");
    save_json(metaJson, meta);

    if args.SaveRolloutLevelMat
        save(rolloutMat, "rollLevel", "-v7.3");
    end

    % ------------------------------------------------------------
    % Terminal summary (quick inspection)
    % ------------------------------------------------------------
    fprintf("\n[A13] Quick pooled summaries (mean across participants)\n");

    % helper: pooled means of participant means (already stored as per-participant means)
    Aglob_MAE = mean(T.A_simple_vs_cpl_global_MAE_mean, "omitnan");
    Apers_MAE = mean(T.A_simple_vs_cpl_person_MAE_mean, "omitnan");
    B_MAE     = mean(T.B_cpl_person_vs_global_MAE_mean, "omitnan");

    Aglob_RMSE = mean(T.A_simple_vs_cpl_global_RMSE_mean, "omitnan");
    Apers_RMSE = mean(T.A_simple_vs_cpl_person_RMSE_mean, "omitnan");
    B_RMSE     = mean(T.B_cpl_person_vs_global_RMSE_mean, "omitnan");

    Aglob_Max = mean(T.A_simple_vs_cpl_global_MaxAbs_mean, "omitnan");
    Apers_Max = mean(T.A_simple_vs_cpl_person_MaxAbs_mean, "omitnan");
    B_Max     = mean(T.B_cpl_person_vs_global_MaxAbs_mean, "omitnan");

    fprintf("  A (simple vs coupled GLOBAL):   MAE_t=%.4g  RMSE_t=%.4g  Max|e|=%.4g\n", Aglob_MAE, Aglob_RMSE, Aglob_Max);
    fprintf("  A (simple vs coupled PERSONAL): MAE_t=%.4g  RMSE_t=%.4g  Max|e|=%.4g\n", Apers_MAE, Apers_RMSE, Apers_Max);
    fprintf("  B (PERSONAL - GLOBAL coupled):  MAE_t=%.4g  RMSE_t=%.4g  Max|e|=%.4g\n", B_MAE, B_RMSE, B_Max);

    % “Is personalization worth it for trust dynamics realism?”
    % Compare divergence to simple: Apers vs Aglob
    if isfinite(Aglob_MAE) && Aglob_MAE > 0 && isfinite(Apers_MAE)
        pct = 100*(Apers_MAE - Aglob_MAE)/Aglob_MAE;
        fprintf("  Personalization effect on A (MAE_t): %+0.2f%% relative to global\n", pct);
    end

    % Fallback audit
    nFallback = sum(T.fallback_applied);
    fprintf("  Guardrail fallback applied for %d/%d participants\n", nFallback, nP);
    if nFallback > 0
        % show brief list
        idx = find(T.fallback_applied);
        for k = 1:min(5,numel(idx))
            i = idx(k);
            fprintf("    pid=%s  used=%s  reason=%s\n", T.participant_id(i), T.model_name_used(i), T.fallback_reason(i));
        end
        if numel(idx) > 5
            fprintf("    ... (%d more)\n", numel(idx)-5);
        end
    end

    fprintf("\n[Step A13] Done.\n");
    fprintf("  Output dir: %s\n", outDir);
    fprintf("  Wrote: %s\n", outCsv);
end

% ======================================================================
% Metrics from an error curve e(t) on the grid
% ======================================================================
function [mae_t, rmse_t, iad, maxabs, acc_abs_curve] = divergence_metrics_from_error(e, dt, Ttotal)
    e = double(e(:));
    e(~isfinite(e)) = 0;

    abs_e = abs(e);

    % integral approximations
    iad = sum(abs_e) * dt;

    mae_t = iad / max(1e-12, Ttotal);
    rmse_t = sqrt( (sum(e.^2) * dt) / max(1e-12, Ttotal) );

    maxabs = max(abs_e);

    % accumulated absolute error curve (trust*s), same length as e
    acc_abs_curve = cumsum(abs_e) * dt;
end

% ======================================================================
% Participant-level summary writing
% ======================================================================
function T = local_write_summary_row(T, rowIdx, X, prefix, qlo, qhi)
    % X: [R x 4] [MAE RMSE IAD MaxAbs]
    prefix = string(prefix);

    mae  = X(:,1); rmse = X(:,2); iad = X(:,3); mx = X(:,4);

    T.(prefix + "_MAE_mean")(rowIdx)  = mean(mae, "omitnan");
    T.(prefix + "_MAE_qlo")(rowIdx)   = quantile_safe(mae, qlo);
    T.(prefix + "_MAE_qhi")(rowIdx)   = quantile_safe(mae, qhi);

    T.(prefix + "_RMSE_mean")(rowIdx) = mean(rmse, "omitnan");
    T.(prefix + "_RMSE_qlo")(rowIdx)  = quantile_safe(rmse, qlo);
    T.(prefix + "_RMSE_qhi")(rowIdx)  = quantile_safe(rmse, qhi);

    T.(prefix + "_IAD_mean")(rowIdx)  = mean(iad, "omitnan");
    T.(prefix + "_IAD_qlo")(rowIdx)   = quantile_safe(iad, qlo);
    T.(prefix + "_IAD_qhi")(rowIdx)   = quantile_safe(iad, qhi);

    T.(prefix + "_MaxAbs_mean")(rowIdx)= mean(mx, "omitnan");
    T.(prefix + "_MaxAbs_qlo")(rowIdx) = quantile_safe(mx, qlo);
    T.(prefix + "_MaxAbs_qhi")(rowIdx) = quantile_safe(mx, qhi);
end

function q = quantile_safe(x, qq)
    x = x(isfinite(x));
    if isempty(x), q = NaN; else, q = quantile(x, qq); end
end

% ======================================================================
% Curve summaries
% ======================================================================
function S = summarize_curves(E_byP, Acc_byP, qlo, qhi)
    % Each is [Ns x nP]
    S = struct();
    S.e_mean = mean(E_byP, 2, "omitnan");
    S.e_qlo  = quantile_cols(E_byP, qlo);
    S.e_qhi  = quantile_cols(E_byP, qhi);

    S.acc_mean = mean(Acc_byP, 2, "omitnan");
    S.acc_qlo  = quantile_cols(Acc_byP, qlo);
    S.acc_qhi  = quantile_cols(Acc_byP, qhi);
end

function q = quantile_cols(X, qq)
    q = NaN(size(X,1),1);
    for i = 1:size(X,1)
        xi = X(i,:);
        xi = xi(isfinite(xi));
        if isempty(xi), q(i) = NaN; else, q(i) = quantile(xi, qq); end
    end
end

% ======================================================================
% Figures
% ======================================================================
function make_fig_error_and_accum(pathPng, pid, t, e, acc, titleStr)
    f = figure('Visible','off');
    % two panels via subplot is acceptable here (single figure file)
    subplot(2,1,1);
    plot(t, e, 'LineWidth', 1.25);
    grid on;
    xlabel('time (s)');
    ylabel('e(t)');
    title(sprintf("pid=%s  %s", string(pid), titleStr), 'Interpreter','none');

    subplot(2,1,2);
    plot(t, acc, 'LineWidth', 1.25);
    grid on;
    xlabel('time (s)');
    ylabel('\int |e| dt (approx)');
    saveas(f, pathPng);
    close(f);
end

function make_fig_pooled_error(pathPng, sgrid, S, titleStr)
    f = figure('Visible','off');
    plot(sgrid, S.e_mean, 'LineWidth', 1.5);
    hold on; grid on;
    plot(sgrid, S.e_qlo, '--', 'LineWidth', 1.0);
    plot(sgrid, S.e_qhi, '--', 'LineWidth', 1.0);
    xlabel('normalized time s \in [0,1]');
    ylabel('e(s)');
    title(titleStr, 'Interpreter','none');
    legend({'mean','qlo','qhi'}, 'Location','best');
    saveas(f, pathPng);
    close(f);
end

function make_fig_pooled_accum(pathPng, sgrid, S, titleStr)
    f = figure('Visible','off');
    plot(sgrid, S.acc_mean, 'LineWidth', 1.5);
    hold on; grid on;
    plot(sgrid, S.acc_qlo, '--', 'LineWidth', 1.0);
    plot(sgrid, S.acc_qhi, '--', 'LineWidth', 1.0);
    xlabel('normalized time s \in [0,1]');
    ylabel('accumulated \int |e| dt (trust*s)');
    title(titleStr, 'Interpreter','none');
    legend({'mean','qlo','qhi'}, 'Location','best');
    saveas(f, pathPng);
    close(f);
end

% ======================================================================
% Global behavior params (reuse A9 logic, but local)
% ======================================================================
function [name, bpar] = resolve_behavior_params_global(modelIdx, fit)
    bpar = struct();
    bpar.tau_flag = 0;
    bpar.m1_flag  = 0;
    bpar.m2_flag  = 0;

    switch modelIdx
        case 0
            name = "model0_trust_as_probability";
            bpar.tau_flag = 1;

        case 1
            name = "model1_threshold";
            assert(isfield(fit,"model1") && isfield(fit.model1,"k_hat"), "[A13] fit.model1.k_hat missing.");
            bpar.m1_flag = 1;
            bpar.k_m1 = fit.model1.k_hat;

        case 2
            name = "model2_offset_lapse";
            assert(isfield(fit,"model2"), "[A13] fit.model2 missing.");
            bpar.m2_flag = 1;
            bpar.k_m2 = fit.model2.k_hat;
            bpar.beta = fit.model2.beta_hat;
            bpar.eps  = fit.model2.eps_hat;

        otherwise
            error("[A13] Unknown model idx: %d", modelIdx);
    end
end

% ======================================================================
% Personalized params resolver (copied from A12; keep consistent)
% ======================================================================
function [bpar, info] = resolve_personalized_behavior_params(pid, Ta10, Tblocks, useA11Guards, epsMax, fallbackStrategy)
    pid = string(pid);
    fallbackStrategy = string(fallbackStrategy);

    idx = find(string(Ta10.participant_id)==pid, 1);
    if isempty(idx)
        error("[A13] pid=%s not found in A10 table.", pid);
    end

    modelA10 = double(Ta10.best_model_idx(idx));
    info = struct();
    info.model_idx_a10 = modelA10;
    info.model_name_a10 = model_name_from_idx(modelA10);

    % Extract parameter estimates from A10 table (support several possible column names)
    kA10    = get_first_existing_numeric(Ta10, idx, ["k_hat","k_mle","k","best_k","k_best"]);
    betaA10 = get_first_existing_numeric(Ta10, idx, ["beta_hat","beta_mle","beta","best_beta","beta_best"]);
    epsA10  = get_first_existing_numeric(Ta10, idx, ["eps_hat","eps_mle","eps","epsilon_hat","epsilon","best_eps","eps_best"]);

    % Under-identified flag from A11 (if available)
    underIdent = false;
    if useA11Guards && ~isempty(Tblocks) && ismember("participant_id", string(Tblocks.Properties.VariableNames))
        j = find(string(Tblocks.participant_id)==pid, 1);
        if ~isempty(j)
            if ismember("flag_under_identified", string(Tblocks.Properties.VariableNames))
                underIdent = logical(Tblocks.flag_under_identified(j));
            elseif ismember("drift_label", string(Tblocks.Properties.VariableNames))
                underIdent = (string(Tblocks.drift_label(j))=="under_identified");
            end
        end
    end

    % Basic parameter sanity checks
    bad2 = ~isfinite(kA10) || kA10 < 0 || ~isfinite(betaA10) || ~isfinite(epsA10) || epsA10 < 0 || epsA10 > epsMax;
    bad1 = ~isfinite(kA10) || kA10 < 0;

    % Decide model_used with fallback
    modelUsed = modelA10;
    fallback_applied = false;
    reason = "";

    if fallbackStrategy ~= "simple"
        warning("[A13] Unknown FallbackStrategy=%s. Using 'simple'.", fallbackStrategy);
    end

    if useA11Guards && underIdent
        fallback_applied = true;
        reason = "A11_under_identified";
        if modelUsed == 2, modelUsed = 1; elseif modelUsed == 1, modelUsed = 0; end
    end

    if modelUsed == 2 && bad2
        fallback_applied = true;
        if strlength(reason)==0, reason="A10_params_invalid_model2"; else, reason=reason + ";A10_params_invalid_model2"; end
        modelUsed = 1;
    end
    if modelUsed == 1 && bad1
        fallback_applied = true;
        if strlength(reason)==0, reason="A10_params_invalid_model1"; else, reason=reason + ";A10_params_invalid_model1"; end
        modelUsed = 0;
    end

    % Build params struct
    bpar = struct();
    bpar.tau_flag = 0;
    bpar.m1_flag  = 0;
    bpar.m2_flag  = 0;

    switch modelUsed
        case 0
            bpar.tau_flag = 1;

        case 1
            bpar.m1_flag = 1;
            bpar.k_m1 = kA10;

        case 2
            bpar.m2_flag = 1;
            bpar.k_m2 = kA10;
            bpar.beta = betaA10;
            bpar.eps  = epsA10;

        otherwise
            error("[A13] Unknown model idx: %d", modelUsed);
    end

    info.model_idx_used = modelUsed;
    info.model_name_used = model_name_from_idx(modelUsed);
    info.fallback_applied = fallback_applied;
    info.fallback_reason = string(reason);
end

function v = get_first_existing_numeric(T, rowIdx, candidates)
    v = NaN;
    candidates = string(candidates);
    for c = candidates
        if ismember(c, string(T.Properties.VariableNames))
            x = T.(c);
            if isnumeric(x)
                v = double(x(rowIdx));
                return;
            end
        end
    end
end

% ======================================================================
% Load theta_star, dt, and VALID participants exactly like A9/A12
% ======================================================================
function [theta_star, dt, participants_valid] = local_load_theta_dt_and_valid_participants_like_A5(run_id)
    run_id = string(run_id);

    % --- A1 archived inputs (VALID participants) ---
    a1Dir = fullfile("derived", "analysis_runs", run_id, "stepA1_prepare_analysis");
    manifestPath = fullfile(a1Dir, "manifest.mat");
    must_exist_file(manifestPath, "A1 manifest");

    validPath = fullfile(a1Dir, "participants_valid_probes_mapped_stepM4.mat");
    must_exist_file(validPath, "A1 VALID participants (mapped probes)");
    participants_valid = load_participants_struct(validPath);

    % --- A3 selection -> resultsMatPath and theta_star ---
    selPath = fullfile("derived","analysis_runs",run_id,"stepA3_model_selection","selection.mat");
    must_exist_file(selPath, "A3 selection.mat");

    Ssel = load(selPath, "selection");
    if ~isfield(Ssel,"selection") || ~isstruct(Ssel.selection)
        error("[A13] A3 selection.mat missing variable 'selection'.");
    end
    selection = Ssel.selection;

    % theta_star
    theta_star = [];
    if isfield(selection,"theta_star") && ~isempty(selection.theta_star)
        theta_star = selection.theta_star(:);
    else
        thetaPath = fullfile("derived","analysis_runs",run_id,"stepA3_model_selection","theta_star.mat");
        if isfile(thetaPath)
            Sth = load(thetaPath);
            theta_star = find_theta_in_struct(Sth);
        end
    end
    if isempty(theta_star)
        error("[A13] Could not resolve theta_star from A3 selection.");
    end

    % results file -> cfg.dt
    if ~isfield(selection,"results_file") || isempty(selection.results_file)
        error("[A13] selection.results_file missing. Cannot locate cfg.dt.");
    end
    resultsMatPath = string(selection.results_file);
    must_exist_file(resultsMatPath, "Fit results MAT (selection.results_file)");

    R = load(resultsMatPath);
    if ~isfield(R,"cfg") || ~isstruct(R.cfg) || ~isfield(R.cfg,"dt") || isempty(R.cfg.dt)
        error("[A13] results MAT does not contain cfg.dt: %s", resultsMatPath);
    end
    dt = double(R.cfg.dt);
    if ~isscalar(dt) || ~isfinite(dt) || dt <= 0
        error("[A13] cfg.dt invalid in results MAT: %s", resultsMatPath);
    end
end

% ======================================================================
% Helpers: participants / tables / alignment / rollout store
% ======================================================================
function ids = local_get_participant_ids(validParticipants)
    if isa(validParticipants, "containers.Map")
        ks = validParticipants.keys;
        ids = string(ks(:));
        ids = sort(ids);
        return;
    end
    if isstruct(validParticipants)
        if ~isfield(validParticipants,"participant_id")
            error("[A13] validParticipants struct missing participant_id field.");
        end
        ids = string({validParticipants.participant_id})';
        ids = sort(ids);
        return;
    end
    if istable(validParticipants)
        if ~ismember("participant_id", string(validParticipants.Properties.VariableNames))
            error("[A13] validParticipants table missing participant_id.");
        end
        ids = string(validParticipants.participant_id);
        ids = sort(ids);
        return;
    end
    error("[A13] Unsupported validParticipants type: %s", class(validParticipants));
end

function T = local_find_first_table(S, preferNames)
    T = [];
    preferNames = string(preferNames);
    for nm = preferNames
        if isfield(S, nm) && istable(S.(nm))
            T = S.(nm);
            return;
        end
    end
    fn = fieldnames(S);
    for k = 1:numel(fn)
        if istable(S.(fn{k}))
            T = S.(fn{k});
            return;
        end
    end
end

function x = local_align_length(x, ref)
    x = double(x(:));
    n = numel(ref);
    if numel(x) == n
        return;
    end
    if numel(x) > n
        x = x(1:n);
    else
        % pad with last value
        if isempty(x)
            x = ref; % worst fallback
        else
            x(end+1:n,1) = x(end);
        end
    end
end

function M = init_rollout_metric_store(nP, R)
    M = struct();
    M.MAE = NaN(nP,R);
    M.RMSE = NaN(nP,R);
    M.IAD = NaN(nP,R);
    M.MaxAbs = NaN(nP,R);
end

% ======================================================================
% Participant collection access (same as A9/A12)
% ======================================================================
function Pp = get_participant_from_collection(collection, pid)
    pid = string(pid);

    if isa(collection, "containers.Map")
        if ~isKey(collection, char(pid))
            error("[A13] VALID participant '%s' not found in collection Map.", pid);
        end
        Pp = collection(char(pid));
        return;
    end

    if isstruct(collection)
        if ~isfield(collection, "participant_id")
            error("[A13] validParticipants struct must have field participant_id.");
        end
        ids = string({collection.participant_id});
        idx = find(ids==pid, 1);
        if isempty(idx), error("[A13] VALID participant '%s' not found in struct collection.", pid); end
        Pp = collection(idx);
        return;
    end

    if istable(collection)
        if ~ismember("participant_id", string(collection.Properties.VariableNames))
            error("[A13] validParticipants table must have participant_id column.");
        end
        idx = find(string(collection.participant_id)==pid, 1);
        if isempty(idx), error("[A13] VALID participant '%s' not found in table collection.", pid); end
        Pp = collection(idx,:);
        return;
    end

    error("[A13] Unsupported validParticipants type: %s", class(collection));
end

function s = sanitize_pid(pid)
    pid = char(string(pid));
    s = regexprep(pid, '[^a-zA-Z0-9_-]', '_');
end

% ======================================================================
% Model name mapping
% ======================================================================
function name = model_name_from_idx(modelIdx)
    switch modelIdx
        case 0, name = "model0_trust_as_probability";
        case 1, name = "model1_threshold";
        case 2, name = "model2_offset_lapse";
        otherwise, name = "unknown";
    end
end
