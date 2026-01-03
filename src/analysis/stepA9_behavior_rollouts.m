function stepA9_behavior_rollouts(run_id, varargin)
% stepA9_behavior_rollouts  Coupled generative rollout analysis on VALID.
%
% Runs CLOSED-LOOP ("coupled") simulations on VALID participants, where door
% decisions are SAMPLED from the behavioral model and then fed back into
% the trust dynamics via the counterfactual outcome inversion rule (handled
% inside trust_simulate_or_predict_one_participant):
%   if sampled_follow ~= recorded_follow  => outcome := 1 - recorded_outcome
%   else                                   outcome := recorded_outcome
%
% This step is NOT about pointwise prediction (A8). It is about whether the
% *coupled* trust+behavior system reproduces interaction *signatures* such
% as override timing and switching structure.
%
% Models (must match A8):
%   Model 0: p_follow = clamp(tau,0,1)                             (no fit)
%   Model 1: p_follow = sigmoid( k*(tau - sc) )                    (k from A8 TRAIN)
%   Model 2: p_follow = (1-eps)*sigmoid(k*tau + beta*(sc-0.5)) + eps*0.5
%
% Outputs written to:
%   derived/analysis_runs/<run_id>/stepA9_behavior_rollouts/
%
% Artifacts:
%   - A9_rollout_stats.mat   (full bundle incl. per-rollout signature arrays)
%   - A9_rollout_stats.csv   (participant-level + pooled summaries)
%   - figures/*.png
%   - meta.mat / meta.json
%
% Signatures (override-centric; robust to high follow base-rate):
%   - follow_rate, override_rate
%   - follow/override rate by block (1..3)   [if block_index available]
%   - switch probabilities: P(switch), P(follow->override), P(override->follow)
%   - inter-override gaps (mean/median/p90)
%   - override streak lengths (mean/p90)
%   - follow streak lengths (mean/p90)  [reported but interpret cautiously]
%
% Requirements on your path:
%   - trust_simulate_or_predict_one_participant (must expose sim.coupled.* in coupled mode)
%   - behavioral_model
%   - must_exist_file, ensure_dir, save_json
%   - load_participants_struct (used in Step A5; preferred). If absent, we fall back to load().
%
% Note:
%   This version DOES NOT duplicate the coupled simulation loop. It reuses
%   trust_simulate_or_predict_one_participant("coupled", ...) and reads
%   sim.coupled.followed_sampled / sim.coupled.block_index / sim.coupled.door_index.

    if nargin < 1 || isempty(run_id)
        error("stepA9_behavior_rollouts: run_id is required.");
    end
    run_id = string(run_id);

    p = inputParser;
    p.addParameter("OutDir", "", @(s) isstring(s) || ischar(s));
    p.addParameter("Overwrite", false, @(x) islogical(x) && isscalar(x));
    p.addParameter("RolloutsPerParticipant", 300, @(x) isnumeric(x) && isscalar(x) && x>=1);
    p.addParameter("Quantiles", [0.05 0.95], @(x) isnumeric(x) && numel(x)==2 && all(x>=0) && all(x<=1));
    p.addParameter("Models", [0 1 2], @(x) isnumeric(x) && isvector(x));
    p.addParameter("RandomSeed", 1, @(x) isnumeric(x) && isscalar(x));
    p.parse(varargin{:});
    args = p.Results;

    rng(args.RandomSeed);

    % ------------------------------------------------------------
    % Load A7 VALID (observed door-level behavior)
    % ------------------------------------------------------------
    a7Dir = fullfile("derived","analysis_runs",run_id,"stepA7_behavior_dataset");
    validMatA7 = fullfile(a7Dir, "behavior_dataset_valid.mat");
    must_exist_file(validMatA7, "A7 VALID dataset");

    S_va = load(validMatA7, "T");
    if ~isfield(S_va,"T") || ~istable(S_va.T), error("[A9] VALID mat missing table T."); end
    Tva = S_va.T;
    Tva = Tva(Tva.is_valid_label==1, :);

    reqCols = ["participant_id","followed","door_index"];
    assert(all(ismember(reqCols, string(Tva.Properties.VariableNames))), "[A9] A7 VALID missing required columns.");
    haveBlockA7 = ismember("block_index", string(Tva.Properties.VariableNames));

    % ------------------------------------------------------------
    % Load A8 fit (behavioral params learned on TRAIN)
    % ------------------------------------------------------------
    a8Dir = fullfile("derived","analysis_runs",run_id,"stepA8_behavior_fit_eval");
    fitMat = fullfile(a8Dir, "fit_params.mat");
    must_exist_file(fitMat, "A8 fit_params.mat");

    S_fit = load(fitMat, "fit");
    if ~isfield(S_fit,"fit"), error("[A9] fit_params.mat missing 'fit' struct."); end
    fit = S_fit.fit;

    % ------------------------------------------------------------
    % Load run-local trust inputs the SAME way Step A5 does:
    %   - VALID participants from A1 archive
    %   - dt from results file referenced by A3 selection
    %   - theta_star from A3 selection
    % ------------------------------------------------------------
    [theta_star, dt, validParticipants] = local_load_theta_dt_and_valid_participants_like_A5(run_id);

    % ------------------------------------------------------------
    % Output directory
    % ------------------------------------------------------------
    outDir = string(args.OutDir);
    if strlength(outDir)==0
        outDir = fullfile("derived","analysis_runs",run_id,"stepA9_behavior_rollouts");
    end
    ensure_dir(outDir);
    figDir = fullfile(outDir, "figures");
    ensure_dir(figDir);

    statsMat = fullfile(outDir, "A9_rollout_stats.mat");
    statsCsv = fullfile(outDir, "A9_rollout_stats.csv");
    metaMat  = fullfile(outDir, "meta.mat");
    metaJson = fullfile(outDir, "meta.json");

    if ~args.Overwrite
        if isfile(statsMat) || isfile(statsCsv)
            error("[A9] Outputs exist. Set Overwrite=true to replace. (%s)", outDir);
        end
    end

    % ------------------------------------------------------------
    % Build observed per-participant sequences (from A7 VALID)
    % ------------------------------------------------------------
    pid_va = string(Tva.participant_id);
    uniqP  = unique(pid_va);
    nP     = numel(uniqP);

    obsStats = table();
    obsStats.participant_id = uniqP;

    for i = 1:nP
        mask = (pid_va == uniqP(i));
        Tpi = Tva(mask, :);
        [~,ord] = sort(double(Tpi.door_index(:)));
        Tpi = Tpi(ord, :);

        seqFollow = double(Tpi.followed(:));      % 1=follow, 0=override
        if haveBlockA7
            seqBlock  = double(Tpi.block_index(:));
        else
            seqBlock  = NaN(size(seqFollow));
        end

        st = compute_signatures(seqFollow, seqBlock);
        obsStats = set_row_struct(obsStats, i, "obs_", st);
    end

    % ------------------------------------------------------------
    % Resolve requested models -> behavior_params structs
    % ------------------------------------------------------------
    models = unique(round(args.Models(:)'));
    R   = args.RolloutsPerParticipant;
    qlo = args.Quantiles(1);
    qhi = args.Quantiles(2);

    modelNames   = strings(1,numel(models));
    modelBParams = cell(1,numel(models));
    for mi = 1:numel(models)
        [modelNames(mi), modelBParams{mi}] = resolve_behavior_params(models(mi), fit);
    end

    % ------------------------------------------------------------
    % Run coupled rollouts via trust_simulate_or_predict_one_participant
    % ------------------------------------------------------------
    roll = struct();
    roll.meta = struct("run_id",char(run_id),"R",R,"models",models,"modelNames",modelNames, ...
                       "qlo",qlo,"qhi",qhi,"random_seed",args.RandomSeed,"dt",dt);
    roll.obsStats = obsStats;

    perModelSummaries = cell(numel(models),1);

    for mi = 1:numel(models)
        mName = modelNames(mi);
        bpar  = modelBParams{mi};

        fprintf("[A9] %s: %d rollouts × %d participants\n", mName, R, nP);

        sig = init_sig_store(nP, R);

        for pi = 1:nP
            pid = uniqP(pi);
            Pp  = get_participant_from_collection(validParticipants, pid);

            for r = 1:R
                % Ensure stochastic rollouts differ but remain reproducible.
                seed = args.RandomSeed + 100000*mi + 1000*pi + r;
                rng(seed);

                sim = trust_simulate_or_predict_one_participant("coupled", theta_star, Pp, dt, bpar);

                if ~isfield(sim, "coupled") || ~isfield(sim.coupled, "followed_sampled")
                    error("[A9] Simulator did not return sim.coupled.followed_sampled in coupled mode (pid=%s).", pid);
                end

                [seqFollow, seqBlock] = coupled_sim_to_sequences(sim);

                st = compute_signatures(seqFollow, seqBlock);
                sig = store_sig(sig, pi, r, st);
            end
        end

        summ = summarize_sig_store(sig, qlo, qhi);
        summ.model = repmat(mName, height(summ), 1);
        summ.participant_id = uniqP;
        summ = join(summ, obsStats, "Keys","participant_id");

        perModelSummaries{mi} = summ;
        roll.(char(mName)).sig = sig;
        roll.(char(mName)).summary = summ;
    end

    allSumm = vertcat(perModelSummaries{:});
    pooled  = pooled_summary(allSumm, modelNames);

    pooled.participant_id = repmat("<POOLED>", height(pooled), 1);
    pooled = movevars(pooled, "participant_id", "Before", 1);

    outTable = [allSumm; pooled];
    writetable(outTable, statsCsv);
    save(statsMat, "roll", "outTable", "pooled", "-v7.3");

    meta = struct();
    meta.run_id = char(run_id);
    meta.created = char(datetime('now','Format','yyyy-MM-dd HH:mm:ss'));
    meta.valid_participants = nP;
    meta.rollouts_per_participant = R;
    meta.models = models;
    meta.model_names = cellstr(modelNames);
    meta.quantiles = [qlo qhi];
    meta.random_seed = args.RandomSeed;
    meta.dt = dt;
    meta.have_block_index_in_A7 = haveBlockA7;
    save(metaMat, "meta");
    save_json(metaJson, meta);

    % ------------------------------------------------------------
    % Figures (pooled; override-centric)
    % ------------------------------------------------------------
    make_fig_pooled_block_rates(fullfile(figDir, "pooled_rates_by_block.png"), outTable, modelNames);
    make_fig_pooled_switch_rates(fullfile(figDir, "pooled_switch_rates.png"), outTable, modelNames);
    make_fig_pooled_override_gap(fullfile(figDir, "pooled_inter_override_gaps.png"), outTable, modelNames);
    make_fig_pooled_override_streak(fullfile(figDir, "pooled_override_streaks.png"), outTable, modelNames);

    fprintf("[Step A9] Done.\n");
    fprintf("  Output dir: %s\n", outDir);
    fprintf("  Wrote: %s\n", statsCsv);
end

% ======================================================================
% Load theta_star, dt, and VALID participants exactly like Step A5
% ======================================================================
function [theta_star, dt, participants_valid] = local_load_theta_dt_and_valid_participants_like_A5(run_id)
    run_id = string(run_id);

    % --- A1 archived inputs (VALID participants) ---
    a1Dir = fullfile("derived", "analysis_runs", run_id, "stepA1_prepare_analysis");
    manifestPath = fullfile(a1Dir, "manifest.mat");
    must_exist_file(manifestPath, "A1 manifest");

    validPath = fullfile(a1Dir, "participants_valid_probes_mapped_stepM4.mat");
    must_exist_file(validPath, "A1 VALID participants (mapped probes)");

    participants_valid = local_load_participants_struct_flexible(validPath);

    % --- A3 selection -> resultsMatPath and theta_star ---
    selPath = fullfile("derived","analysis_runs",run_id,"stepA3_model_selection","selection.mat");
    must_exist_file(selPath, "A3 selection.mat");

    Ssel = load(selPath, "selection");
    if ~isfield(Ssel,"selection") || ~isstruct(Ssel.selection)
        error("[A9] A3 selection.mat missing variable 'selection'.");
    end
    selection = Ssel.selection;

    % theta_star
    theta_star = [];
    if isfield(selection,"theta_star") && ~isempty(selection.theta_star)
        theta_star = selection.theta_star(:);
    else
        % fallback: A3 theta_star.mat
        thetaPath = fullfile("derived","analysis_runs",run_id,"stepA3_model_selection","theta_star.mat");
        if isfile(thetaPath)
            Sth = load(thetaPath);
            theta_star = local_find_theta_in_struct(Sth);
        end
    end
    if isempty(theta_star)
        error("[A9] Could not resolve theta_star from A3 selection (selection.theta_star missing/empty).");
    end

    % results file -> cfg.dt
    if ~isfield(selection,"results_file") || isempty(selection.results_file)
        error("[A9] selection.results_file missing. Cannot locate cfg.dt.");
    end
    resultsMatPath = string(selection.results_file);
    must_exist_file(resultsMatPath, "Fit results MAT (selection.results_file)");

    R = load(resultsMatPath);
    if ~isfield(R,"cfg") || ~isstruct(R.cfg) || ~isfield(R.cfg,"dt") || isempty(R.cfg.dt)
        error("[A9] results MAT does not contain cfg.dt: %s", resultsMatPath);
    end
    dt = double(R.cfg.dt);
    if ~isscalar(dt) || ~isfinite(dt) || dt <= 0
        error("[A9] cfg.dt invalid in results MAT: %s", resultsMatPath);
    end
end

function P = local_load_participants_struct_flexible(matPath)
    % Prefer project helper (used in A5). Fallback to common variable names.
    if exist("load_participants_struct","file") == 2
        P = load_participants_struct(matPath);
        return;
    end

    S = load(matPath);
    % Common conventions:
    %   - participants (struct array)
    %   - P (struct array)
    %   - participants_valid
    cand = ["participants","P","participants_valid","participants_train"];
    for c = cand
        if isfield(S, c)
            P = S.(c);
            if isstruct(P)
                return;
            end
        end
    end
    error("[A9] Could not load participants struct from %s (no known variable names).", matPath);
end

function theta = local_find_theta_in_struct(S)
    % Minimal version of A5's find_theta_in_struct.
    theta = [];
    if isempty(S) || ~isstruct(S), return; end
    f = fieldnames(S);
    for i = 1:numel(f)
        v = S.(f{i});
        if isnumeric(v) && (isvector(v) || ismatrix(v)) && numel(v) >= 1
            % Heuristic: accept first numeric vector-ish thing
            theta = v(:);
            return;
        end
        if isstruct(v)
            theta = local_find_theta_in_struct(v);
            if ~isempty(theta), return; end
        end
    end
end

% ======================================================================
% Convert coupled sim output -> sequences (robust to missing block index)
% ======================================================================
function [seqFollow, seqBlock] = coupled_sim_to_sequences(sim)
    seqFollow = double(sim.coupled.followed_sampled(:));  % 1=follow, 0=override

    % door_index sorting if present
    if isfield(sim.coupled, "door_index") && ~isempty(sim.coupled.door_index) && any(isfinite(double(sim.coupled.door_index(:))))
        doorIdx = double(sim.coupled.door_index(:));
        [~,ord] = sort(doorIdx);
        seqFollow = seqFollow(ord);

        if isfield(sim.coupled, "block_index") && ~isempty(sim.coupled.block_index)
            seqBlock = double(sim.coupled.block_index(:));
            seqBlock = seqBlock(ord);
        else
            seqBlock = NaN(size(seqFollow));
        end
    else
        if isfield(sim.coupled, "block_index") && ~isempty(sim.coupled.block_index)
            seqBlock = double(sim.coupled.block_index(:));
        else
            seqBlock = NaN(size(seqFollow));
        end
    end

    if any(isnan(seqFollow))
        seqFollow(isnan(seqFollow)) = 1; % keep run alive, but clearly degenerate
    end
end

% ======================================================================
% Build behavior_params for your behavioral_model(state, params)
% ======================================================================
function [name, bpar] = resolve_behavior_params(modelIdx, fit)
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
            assert(isfield(fit,"model1") && isfield(fit.model1,"k_hat"), "[A9] fit.model1.k_hat missing.");
            k = fit.model1.k_hat;
            bpar.m1_flag = 1;
            bpar.k_m1 = k;

        case 2
            name = "model2_offset_lapse";
            assert(isfield(fit,"model2"), "[A9] fit.model2 missing.");
            k    = fit.model2.k_hat;
            beta = fit.model2.beta_hat;
            eps  = fit.model2.eps_hat;

            bpar.m2_flag = 1;
            bpar.k_m2 = k;
            bpar.beta = beta;
            bpar.eps  = eps;

        otherwise
            error("[A9] Unknown model index: %d", modelIdx);
    end
end

% ======================================================================
% Signatures (override-centric + transitions)
% ======================================================================
function st = compute_signatures(seqFollow, seqBlock)
    seqFollow = double(seqFollow(:));
    seqOverride = 1 - seqFollow;

    n = numel(seqFollow);
    if n==0
        st = empty_sig();
        return;
    end

    st.N_doors = n;

    st.follow_rate = mean(seqFollow);
    st.override_rate = mean(seqOverride);

    if all(isnan(seqBlock))
        st.follow_rate_b1 = NaN; st.follow_rate_b2 = NaN; st.follow_rate_b3 = NaN;
        st.override_rate_b1 = NaN; st.override_rate_b2 = NaN; st.override_rate_b3 = NaN;
    else
        st.follow_rate_b1 = rate_in_block(seqFollow, seqBlock, 1);
        st.follow_rate_b2 = rate_in_block(seqFollow, seqBlock, 2);
        st.follow_rate_b3 = rate_in_block(seqFollow, seqBlock, 3);

        st.override_rate_b1 = rate_in_block(seqOverride, seqBlock, 1);
        st.override_rate_b2 = rate_in_block(seqOverride, seqBlock, 2);
        st.override_rate_b3 = rate_in_block(seqOverride, seqBlock, 3);
    end

    if n >= 2
        sw = (seqFollow(2:end) ~= seqFollow(1:end-1));
        st.p_switch = mean(sw);

        f2o = (seqFollow(1:end-1)==1) & (seqFollow(2:end)==0);
        o2f = (seqFollow(1:end-1)==0) & (seqFollow(2:end)==1);
        st.p_follow_to_override = mean(f2o);
        st.p_override_to_follow = mean(o2f);
    else
        st.p_switch = NaN;
        st.p_follow_to_override = NaN;
        st.p_override_to_follow = NaN;
    end

    ovIdx = find(seqOverride==1);
    if numel(ovIdx) >= 2
        gaps = diff(ovIdx);
        st.inter_override_gap_mean = mean(gaps);
        st.inter_override_gap_med  = median(gaps);
        st.inter_override_gap_p90  = prctile(gaps, 90);
    else
        st.inter_override_gap_mean = NaN;
        st.inter_override_gap_med  = NaN;
        st.inter_override_gap_p90  = NaN;
    end

    st.override_streak_mean = mean_streak_len(seqOverride);
    st.override_streak_p90  = prctile_streak_len(seqOverride, 90);

    st.follow_streak_mean = mean_streak_len(seqFollow);
    st.follow_streak_p90  = prctile_streak_len(seqFollow, 90);
end

function st = empty_sig()
    st = struct( ...
        "N_doors",0, ...
        "follow_rate",NaN,"override_rate",NaN, ...
        "follow_rate_b1",NaN,"follow_rate_b2",NaN,"follow_rate_b3",NaN, ...
        "override_rate_b1",NaN,"override_rate_b2",NaN,"override_rate_b3",NaN, ...
        "p_switch",NaN,"p_follow_to_override",NaN,"p_override_to_follow",NaN, ...
        "inter_override_gap_mean",NaN,"inter_override_gap_med",NaN,"inter_override_gap_p90",NaN, ...
        "override_streak_mean",NaN,"override_streak_p90",NaN, ...
        "follow_streak_mean",NaN,"follow_streak_p90",NaN ...
    );
end

function r = rate_in_block(x, b, blk)
    mask = (double(b(:)) == blk);
    if any(mask), r = mean(double(x(mask))); else, r = NaN; end
end

function m = mean_streak_len(binarySeq)
    lens = streak_lengths(binarySeq);
    if isempty(lens), m = 0; else, m = mean(lens); end
end

function p = prctile_streak_len(binarySeq, q)
    lens = streak_lengths(binarySeq);
    if isempty(lens), p = 0; else, p = prctile(lens, q); end
end

function lens = streak_lengths(binarySeq)
    x = double(binarySeq(:)~=0);
    if isempty(x), lens = []; return; end
    d = diff([0; x; 0]);
    runStarts = find(d==1);
    runEnds   = find(d==-1) - 1;
    lens = runEnds - runStarts + 1;
end

% ======================================================================
% Storage and summarization
% ======================================================================
function sig = init_sig_store(nP, R)
    fields = fieldnames(empty_sig());
    for i = 1:numel(fields)
        sig.(fields{i}) = NaN(nP, R);
    end
end

function sig = store_sig(sig, pi, r, st)
    fn = fieldnames(st);
    for i = 1:numel(fn)
        sig.(fn{i})(pi, r) = st.(fn{i});
    end
end

function T = summarize_sig_store(sig, qlo, qhi)
    fn = fieldnames(sig);
    T = table();
    for i = 1:numel(fn)
        X = sig.(fn{i}); % [nP x R]
        T.(fn{i} + "_mean") = mean(X, 2, "omitnan");
        T.(fn{i} + "_qlo")  = quantile_rows(X, qlo);
        T.(fn{i} + "_qhi")  = quantile_rows(X, qhi);
    end
end

function q = quantile_rows(X, qq)
    q = NaN(size(X,1),1);
    for i = 1:size(X,1)
        xi = X(i,:);
        xi = xi(isfinite(xi));
        if isempty(xi), q(i) = NaN; else, q(i) = quantile(xi, qq); end
    end
end

function T = set_row_struct(T, rowIdx, prefix, st)
    fn = fieldnames(st);
    for i = 1:numel(fn)
        col = prefix + string(fn{i});
        if ~ismember(col, string(T.Properties.VariableNames))
            T.(col) = NaN(height(T),1);
        end
        T.(col)(rowIdx) = st.(fn{i});
    end
end

function pooled = pooled_summary(allSumm, modelNames)
    % Build a pooled table with the SAME variables (and types) as allSumm,
    % except participant_id (which is handled outside).
    %
    % We take a 1-row template from allSumm, strip participant_id, and then
    % overwrite numeric columns with pooled means per model.

    % Variables to keep in pooled (everything except participant_id)
    varsKeep = allSumm.Properties.VariableNames;
    varsKeep = varsKeep(~strcmp(varsKeep, "participant_id"));

    % Prepare an empty pooled table with identical vars/types to allSumm(varsKeep)
    tmpl = allSumm(1, varsKeep);
    tmpl(1,:) = [];  % make empty but keep schema
    pooled = tmpl;

    % Identify numeric/logical columns to pool (skip "model" and any non-numeric)
    numericVars = {};
    for i = 1:numel(varsKeep)
        v = varsKeep{i};
        if strcmp(v, "model")
            continue;
        end
        x = allSumm.(v);
        if isnumeric(x) || islogical(x)
            numericVars{end+1} = v; %#ok<AGROW>
        end
    end

    % Build one pooled row per model with full schema
    for mi = 1:numel(modelNames)
        m = string(modelNames(mi));

        % Subset rows for this model (exclude <POOLED> if any already exists)
        M = allSumm(string(allSumm.model)==m & string(allSumm.participant_id)~="<POOLED>", :);

        % Start from a 1-row template that already has all columns/types
        row = allSumm(1, varsKeep);

        % Set model
        row.model = m;

        % Fill pooled means for numeric/logical columns
        for ci = 1:numel(numericVars)
            c = numericVars{ci};
            row.(c) = mean(M.(c), "omitnan");
        end

        pooled = [pooled; row]; %#ok<AGROW>
    end
end


% ======================================================================
% Plotting (simple pooled plots)
% ======================================================================
function make_fig_pooled_block_rates(pathPng, outTable, modelNames)
    f = figure('Visible','off');
    x = [1 2 3];
    hold on;

    for mi = 1:numel(modelNames)
        m = modelNames(mi);
        M = outTable(string(outTable.model)==string(m) & outTable.participant_id~="<POOLED>", :);

        y = [ ...
            mean(M.follow_rate_b1_mean,'omitnan'), ...
            mean(M.follow_rate_b2_mean,'omitnan'), ...
            mean(M.follow_rate_b3_mean,'omitnan') ];

        if all(~isfinite(y)), continue; end
        plot(x, y, '-o', 'LineWidth', 1.5, 'DisplayName', char(m));
    end

    Mop = outTable(string(outTable.model)==string(modelNames(1)) & outTable.participant_id~="<POOLED>", :);
    yobs = [ ...
        mean(Mop.obs_follow_rate_b1,'omitnan'), ...
        mean(Mop.obs_follow_rate_b2,'omitnan'), ...
        mean(Mop.obs_follow_rate_b3,'omitnan') ];

    if any(isfinite(yobs))
        plot(x, yobs, '--s', 'LineWidth', 1.5, 'DisplayName', 'observed');
    end

    xlabel('Block');
    ylabel('Follow rate');
    title('Pooled follow rate by block (observed vs coupled rollouts)');
    grid on;
    xlim([1 3]);
    ylim([0 1]);
    legend('Location','best');
    saveas(f, pathPng);
    close(f);
end

function make_fig_pooled_switch_rates(pathPng, outTable, modelNames)
    f = figure('Visible','off');

    cats = ["p_switch_mean","p_follow_to_override_mean","p_override_to_follow_mean"];
    catLabels = {'P(switch)','P(F→O)','P(O→F)'};
    x = 1:numel(cats);

    hold on;
    for mi = 1:numel(modelNames)
        m = modelNames(mi);
        M = outTable(string(outTable.model)==string(m) & outTable.participant_id~="<POOLED>", :);
        y = zeros(1,numel(cats));
        for ci = 1:numel(cats)
            y(ci) = mean(M.(cats(ci)),'omitnan');
        end
        plot(x, y, '-o', 'LineWidth', 1.5, 'DisplayName', char(m));
    end

    Mop = outTable(string(outTable.model)==string(modelNames(1)) & outTable.participant_id~="<POOLED>", :);
    yobs = [ ...
        mean(Mop.obs_p_switch,'omitnan'), ...
        mean(Mop.obs_p_follow_to_override,'omitnan'), ...
        mean(Mop.obs_p_override_to_follow,'omitnan') ];
    plot(x, yobs, '--s', 'LineWidth', 1.5, 'DisplayName', 'observed');

    set(gca,'XTick',x,'XTickLabel',catLabels);
    ylabel('Probability');
    title('Pooled transition structure (observed vs coupled rollouts)');
    grid on;
    ylim([0 1]);
    legend('Location','best');
    saveas(f, pathPng);
    close(f);
end

function make_fig_pooled_override_gap(pathPng, outTable, modelNames)
    f = figure('Visible','off');
    hold on;

    for mi = 1:numel(modelNames)
        m = modelNames(mi);
        M = outTable(string(outTable.model)==string(m) & outTable.participant_id~="<POOLED>", :);
        y = mean(M.inter_override_gap_mean_mean,'omitnan');
        plot(mi, y, 'o', 'MarkerSize', 8, 'DisplayName', char(m));
    end

    Mop = outTable(string(outTable.model)==string(modelNames(1)) & outTable.participant_id~="<POOLED>", :);
    yobs = mean(Mop.obs_inter_override_gap_mean,'omitnan');
    plot(numel(modelNames)+1, yobs, 's', 'MarkerSize', 8, 'DisplayName', 'observed');

    set(gca,'XTick',1:(numel(modelNames)+1), 'XTickLabel',[cellstr(modelNames) {'observed'}]);
    ylabel('Mean gap (doors) between overrides');
    title('Pooled inter-override gap (override sparsity / clustering)');
    grid on;
    legend('Location','best');
    saveas(f, pathPng);
    close(f);
end

function make_fig_pooled_override_streak(pathPng, outTable, modelNames)
    f = figure('Visible','off');
    hold on;

    for mi = 1:numel(modelNames)
        m = modelNames(mi);
        M = outTable(string(outTable.model)==string(m) & outTable.participant_id~="<POOLED>", :);
        y = mean(M.override_streak_mean_mean,'omitnan');
        plot(mi, y, 'o', 'MarkerSize', 8, 'DisplayName', char(m));
    end

    Mop = outTable(string(outTable.model)==string(modelNames(1)) & outTable.participant_id~="<POOLED>", :);
    yobs = mean(Mop.obs_override_streak_mean,'omitnan');
    plot(numel(modelNames)+1, yobs, 's', 'MarkerSize', 8, 'DisplayName', 'observed');

    set(gca,'XTick',1:(numel(modelNames)+1), 'XTickLabel',[cellstr(modelNames) {'observed'}]);
    ylabel('Mean override streak length');
    title('Pooled override streak length (burstiness)');
    grid on;
    legend('Location','best');
    saveas(f, pathPng);
    close(f);
end

% ======================================================================
% Participant collection access (struct / table / map)
% ======================================================================
function Pp = get_participant_from_collection(collection, pid)
    if isa(collection, "containers.Map")
        if ~isKey(collection, char(pid))
            error("[A9] VALID participant '%s' not found in collection Map.", pid);
        end
        Pp = collection(char(pid));
        return;
    end

    if isstruct(collection)
        if ~isfield(collection, "participant_id")
            error("[A9] validParticipants struct must have field participant_id.");
        end
        ids = string({collection.participant_id});
        idx = find(ids==pid, 1);
        if isempty(idx), error("[A9] VALID participant '%s' not found in struct collection.", pid); end
        Pp = collection(idx);
        return;
    end

    if istable(collection)
        if ~ismember("participant_id", string(collection.Properties.VariableNames))
            error("[A9] validParticipants table must have participant_id column.");
        end
        idx = find(string(collection.participant_id)==pid, 1);
        if isempty(idx), error("[A9] VALID participant '%s' not found in table collection.", pid); end
        Pp = collection(idx,:);
        return;
    end

    error("[A9] Unsupported validParticipants type: %s", class(collection));
end
