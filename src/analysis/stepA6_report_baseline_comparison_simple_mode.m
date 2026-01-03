function stepA6_report_baseline_comparison_simple_mode(run_id, inA5Dir)
% stepA6_report_baseline_comparison_simple_mode  Generate thesis-ready reports and figures from Step A5 outputs.
%
% This step consumes the run-local plot bundle produced by Step A5
% (A5_plot_data.mat) and generates publication-ready tables and figures for the
% baseline comparison in SIMPLE mode. A6 intentionally avoids recomputing
% metrics from CSV files; instead it loads the consolidated MAT bundle to ensure
% consistent reporting and reduce numerical drift.
%
% In addition to aggregate tables/figures, this step generates per-participant
% trajectory plots on the validation set using the same simulation helpers that
% produced the aligned measurements in A5:
%   - trust_simulate_or_predict_one_participant("simple", ...)
%   - trust_simulate_baseline_one_participant(method, ...)
%
% Inputs
%   run_id (string/char)
%       Analysis run identifier. Outputs are written under:
%         derived/analysis_runs/<run_id>/stepA6_report/
%
%   inA5Dir (string/char) (optional)
%       Folder containing Step A5 outputs. If omitted or empty, defaults to:
%         derived/analysis_runs/<run_id>/stepA5_baseline_comparison
%
% Outputs
%   (none)
%       Files written to:
%         derived/analysis_runs/<run_id>/stepA6_report/
%           meta.mat
%           A6_valid_overall_table.csv
%           A6_valid_kind_table.csv
%           A6_valid_improvements_table.csv
%           A6_best_baseline.txt
%           fig_valid_overall_wrmse.png
%           fig_valid_wrmse_by_kind.png
%           fig_participant_wrmse_model_vs_bestbaseline.png   (best-effort)
%           participant_trajectories/
%
% UPDATED (trajectory organization + diagnostics + streamlined figures):
%   - Plot trajectories for ALL validation participants (no cherry picking).
%   - Store plots under participant_trajectories/<participant_id>/ subfolders.
%   - Per participant: generate ONLY TWO figures:
%       (1) A 3x2 multi-panel grid with FIVE methods (one subplot each),
%           with the bottom-right (6th) panel reserved for a wRMSE note box.
%       (2) A combo plot (constants dashed + model + best baseline) kept as-is.
%
% Assumptions / Dependencies
%   - Step A5 has produced <A5Dir>/A5_plot_data.mat containing variable A5plot.
%   - Step A1 archived validation participants and weights exist at:
%       derived/analysis_runs/<run_id>/stepA1_prepare_analysis/
%         participants_valid_probes_mapped_stepM4.mat
%         measurement_weights.mat
%   - Utility functions are available on the MATLAB path:
%       ensure_dir, load_participants_struct, read_participant_ids
%   - Simulation helper functions are available on the MATLAB path:
%       trust_simulate_or_predict_one_participant, trust_simulate_baseline_one_participant
%
% Notes
%   - Baseline naming: method "bump_asymmetric" refers to the literature-style
%     bump+saturation baseline in discrete-time form. The method string is kept
%     for compatibility with Step A5 artifacts.
%   - Plotting: trajectory plots include door event markers and measurement
%     markers overlaid on model/baseline trust trajectories.

    if nargin < 1 || isempty(run_id)
        error("stepA6_report_baseline_comparison_simple_mode: run_id is required.");
    end
    run_id = string(run_id);

    if nargin < 2 || isempty(inA5Dir)
        inA5Dir = fullfile("derived","analysis_runs",run_id,"stepA5_baseline_comparison");
    end
    inA5Dir = string(inA5Dir);

    if ~isfolder(inA5Dir)
        error("A5 folder not found: %s", inA5Dir);
    end

    % ------------------------------------------------------------------
    % Load consolidated A5 plot bundle (single source of truth for tables)
    % ------------------------------------------------------------------
    plotMat = fullfile(inA5Dir, "A5_plot_data.mat");
    if ~isfile(plotMat)
        error("A5_plot_data.mat not found. Re-run updated Step A5. Missing: %s", plotMat);
    end

    S = load(plotMat, "A5plot");
    if ~isfield(S, "A5plot")
        error("A5_plot_data.mat does not contain variable 'A5plot'.");
    end
    A5plot = S.A5plot;

    % ------------------------------------------------------------------
    % Output directory
    % ------------------------------------------------------------------
    outDir = fullfile("derived","analysis_runs",run_id,"stepA6_report");
    if exist("ensure_dir","file") == 2
        ensure_dir(outDir);
    else
        if ~isfolder(outDir), mkdir(outDir); end
    end

    % Pull core tables/bundles from A5 plot data
    methods      = string(A5plot.methods);
    validOverall = A5plot.validOverall;
    validKind    = A5plot.validKind;
    improvements = A5plot.improvements;

    bestBaseline = string(A5plot.bestBaseline_valid);
    bestVal      = double(A5plot.bestBaseline_valid_wRMSE);

    % ---- Sanity checks ----
    if ~ismember(bestBaseline, methods) && ~(ismissing(bestBaseline) || strlength(bestBaseline)==0)
        warning("[A6] bestBaseline '%s' not found in A5plot.methods.", bestBaseline);
    end

    % Reorder validOverall to match methods ordering (if needed)
    if istable(validOverall) && ismember("method", validOverall.Properties.VariableNames)
        vMethods = string(validOverall.method);
        if ~isequal(vMethods(:), methods(:))
            [tf, loc] = ismember(methods, vMethods);
            if all(tf)
                validOverall = validOverall(loc, :);
            else
                warning("[A6] Could not fully reorder validOverall to match methods.");
            end
        end
    end

    % ------------------------------------------------------------------
    % Export thesis tables as CSV (from MAT bundle; no recomputation)
    % ------------------------------------------------------------------
    writetable(validOverall, fullfile(outDir, "A6_valid_overall_table.csv"));
    writetable(validKind,    fullfile(outDir, "A6_valid_kind_table.csv"));

    if ~isempty(improvements) && istable(improvements) && height(improvements) > 0
        impValid = improvements(string(improvements.split) == "valid", :);
        writetable(impValid, fullfile(outDir, "A6_valid_improvements_table.csv"));
    else
        impValid = table(); %#ok<NASGU>
        warning("[A6] No improvements table available in A5plot bundle.");
    end

    % Record best baseline choice for quick reference (human-readable)
    fid = fopen(fullfile(outDir, "A6_best_baseline.txt"), "w");
    fprintf(fid, "Best baseline on VALID (overall wRMSE): %s (wRMSE=%.6g)\n", bestBaseline, bestVal);
    fclose(fid);

    % ------------------------------------------------------------------
    % Meta information (provenance for reporting)
    % ------------------------------------------------------------------
    meta = struct();
    meta.run_id = char(run_id);
    meta.inA5Dir = char(inA5Dir);
    meta.outDir = char(outDir);
    meta.created = char(datetime('now','Format','yyyy-MM-dd HH:mm:ss'));
    meta.bestBaseline_valid = char(bestBaseline);
    meta.bestBaseline_valid_wRMSE = bestVal;
    meta.source_plot_bundle = char(plotMat);
    meta.note = "A6 uses A5_plot_data.mat for tables; trajectories re-load VALID participants from A1 archive and simulate via trust_simulate_* functions.";
    save(fullfile(outDir, "meta.mat"), "meta");

    % ------------------------------------------------------------------
    % FIG 1: VALID overall wRMSE by method
    % ------------------------------------------------------------------
    fig1 = figure('Visible','off','Color','w','Name','VALID overall wRMSE');
    vals = validOverall.wRMSE_overall;
    cats = categorical(string(validOverall.method), cellstr(methods), cellstr(methods));
    bar(cats, vals);
    grid on;
    ylabel("wRMSE (VALID, overall)");
    title("Validation overall wRMSE by method");
    ax = gca;
    ax.XTickLabelRotation = 25;
    exportgraphics(fig1, fullfile(outDir, "fig_valid_overall_wrmse.png"), 'Resolution', 220);
    close(fig1);

    % ------------------------------------------------------------------
    % FIG 2: VALID wRMSE by measurement kind (grouped bar; union of kinds)
    % ------------------------------------------------------------------
    if ~istable(validKind) || height(validKind) == 0
        warning("[A6] validKind table missing/empty; skipping fig_valid_wrmse_by_kind.png");
    else
        mk = string(validKind.method) + "||" + string(validKind.kind);
        [u, ~, ic] = unique(mk);
        counts = accumarray(ic, 1);
        if any(counts > 1)
            bad = u(counts > 1);
            error("[A6] validKind has duplicate (method,kind) rows. Examples:\n  %s", strjoin(bad(1:min(10,numel(bad))), "\n  "));
        end

        kindsAll = sort(unique(string(validKind.kind)));
        M = numel(methods);
        K = numel(kindsAll);
        matWR = NaN(M, K);

        for i = 1:M
            m = methods(i);
            for k = 1:K
                kk = kindsAll(k);
                row = validKind(string(validKind.method) == m & string(validKind.kind) == kk, :);
                if ~isempty(row)
                    matWR(i,k) = row.wRMSE(1);
                end
            end
        end

        fig2 = figure('Visible','off','Color','w','Name','VALID wRMSE by kind');
        bar(matWR, 'grouped');
        grid on;
        xlabel("method");
        ylabel("wRMSE (VALID)");
        title("Validation wRMSE by measurement type (kind)");
        ax = gca;
        ax.XTick = 1:M;
        ax.XTickLabel = cellstr(methods);
        ax.XTickLabelRotation = 25;
        legend(cellstr(kindsAll), 'Location','bestoutside');
        exportgraphics(fig2, fullfile(outDir, "fig_valid_wrmse_by_kind.png"), 'Resolution', 220);
        close(fig2);
    end

    % ------------------------------------------------------------------
    % FIG 3: Participant-level wRMSE distribution (model vs best baseline)
    % ------------------------------------------------------------------
    Tjoin = table();
    canFig3 = true;

    if ~isfield(A5plot, "allParticipants") || ~isstruct(A5plot.allParticipants), canFig3 = false; end
    if canFig3 && ~isfield(A5plot.allParticipants, "simple_selected"), canFig3 = false; end
    if canFig3 && ~isfield(A5plot.allParticipants.simple_selected, "valid"), canFig3 = false; end
    if canFig3 && ~(isfield(A5plot.allParticipants, bestBaseline) && isfield(A5plot.allParticipants.(bestBaseline), "valid"))
        canFig3 = false;
    end

    if ~canFig3 || ismissing(bestBaseline) || strlength(bestBaseline)==0
        warning("[A6] Missing participant tables for model/bestBaseline; skipping fig_participant_wrmse_model_vs_bestbaseline.png");
    else
        Tp_model = A5plot.allParticipants.simple_selected.valid;
        Tp_base  = A5plot.allParticipants.(bestBaseline).valid;

        if ~istable(Tp_model) || ~istable(Tp_base)
            warning("[A6] Participant tables are not tables; skipping fig 3.");
        else
            mustCols = ["participant_id","wRMSE","N"];
            for c = mustCols
                if ~ismember(c, string(Tp_model.Properties.VariableNames))
                    error("Tp_model missing required column '%s'.", c);
                end
                if ~ismember(c, string(Tp_base.Properties.VariableNames))
                    error("Tp_base missing required column '%s'.", c);
                end
            end

            Tp_model.participant_id = string(Tp_model.participant_id);
            Tp_base.participant_id  = string(Tp_base.participant_id);

            Tp_model = Tp_model(:, {'participant_id','wRMSE','N'});
            Tp_base  = Tp_base(:,  {'participant_id','wRMSE','N'});

            Tp_model = local_rename_var(Tp_model, "wRMSE", "wRMSE_model");
            Tp_base  = local_rename_var(Tp_base,  "wRMSE", "wRMSE_baseline");
            Tp_model = local_rename_var(Tp_model, "N",     "N_model");
            Tp_base  = local_rename_var(Tp_base,  "N",     "N_baseline");

            Tjoin = innerjoin(Tp_model, Tp_base, 'Keys', 'participant_id');

            fig3 = figure('Visible','off','Color','w','Name','Participant wRMSE distribution');
            data = [Tjoin.wRMSE_model, Tjoin.wRMSE_baseline];
            boxplot(data, 'Labels', {char("simple\_selected"), char(bestBaseline)});
            grid on;
            ylabel("participant-level wRMSE (VALID)");
            title("Participant-level validation wRMSE: model vs best baseline");
            exportgraphics(fig3, fullfile(outDir, "fig_participant_wrmse_model_vs_bestbaseline.png"), 'Resolution', 220);
            close(fig3);
        end
    end

    % ------------------------------------------------------------------
    % Participant trajectory plots (validation set; ALL participants)
    %   - 2 figures per participant:
    %       (A) methods grid (5 subplots + note panel)
    %       (B) combo plot (kept as-is)
    % ------------------------------------------------------------------
    trajDir = fullfile(outDir, "participant_trajectories");
    if exist("ensure_dir","file") == 2
        ensure_dir(trajDir);
    else
        if ~isfolder(trajDir), mkdir(trajDir); end
    end

    a1Dir = fullfile("derived","analysis_runs",run_id,"stepA1_prepare_analysis");
    validMat = fullfile(a1Dir, "participants_valid_probes_mapped_stepM4.mat");
    weightsMat = fullfile(a1Dir, "measurement_weights.mat");

    if ~isfield(A5plot, "dt") || isempty(A5plot.dt)
        warning("[A6] A5plot.dt missing; skipping trajectory plots.");
    elseif ~isfield(A5plot, "fit_params") || ~isfield(A5plot.fit_params, "globalConstTrain")
        warning("[A6] A5plot.fit_params missing; skipping trajectory plots.");
    elseif ~isfile(validMat)
        warning("[A6] VALID participants file not found: %s. Skipping trajectory plots.", validMat);
    elseif ~isfile(weightsMat)
        warning("[A6] Weights file not found: %s. Skipping trajectory plots.", weightsMat);
    else
        participants_valid = load_participants_struct(validMat);

        W = load(weightsMat, "weights");
        if ~isfield(W, "weights")
            warning("[A6] weightsMat missing variable 'weights'; skipping trajectory plots.");
        else
            weights = W.weights;
            dt = double(A5plot.dt);

            theta_star = [];
            if isfield(A5plot, "meta") && isfield(A5plot.meta, "theta_star")
                theta_star = A5plot.meta.theta_star(:);
            elseif isfield(A5plot, "theta_star")
                theta_star = A5plot.theta_star(:);
            end
            if isempty(theta_star)
                warning("[A6] theta_star not found in A5plot bundle; skipping trajectory plots.");
            else
                % ---- cfg: carry clip01 + optimo knobs ----
                cfg = struct('clip01', true);
                if isfield(A5plot, "cfg") && isfield(A5plot.cfg, "clip01")
                    cfg.clip01 = logical(A5plot.cfg.clip01);
                elseif isfield(A5plot, "meta") && isfield(A5plot.meta, "cfg") && isfield(A5plot.meta.cfg, "clip01")
                    cfg.clip01 = logical(A5plot.meta.cfg.clip01);
                end

                if isfield(A5plot, "cfg") && isfield(A5plot.cfg, "optimo")
                    cfg.optimo = A5plot.cfg.optimo;
                elseif isfield(A5plot, "meta") && isfield(A5plot.meta, "cfg") && isfield(A5plot.meta.cfg, "optimo")
                    cfg.optimo = A5plot.meta.cfg.optimo;
                end

                % ---- baselineParams: include optimo params too ----
                fit_params = A5plot.fit_params;
                baselineParams = struct();
                baselineParams.globalConstTrain = double(fit_params.globalConstTrain);

                dp = NaN; dm = NaN;
                if isfield(fit_params, "ode") && isstruct(fit_params.ode)
                    if isfield(fit_params.ode, "delta_plus"),  dp = double(fit_params.ode.delta_plus); end
                    if isfield(fit_params.ode, "delta_minus"), dm = double(fit_params.ode.delta_minus); end
                end
                if ~isfinite(dp), dp = 0; end
                if ~isfinite(dm), dm = 0; end
                baselineParams.delta_plus  = max(dp, 0);
                baselineParams.delta_minus = max(dm, 0);

                if isfield(fit_params, "bump") && isstruct(fit_params.bump) && isfield(fit_params.bump, "delta")
                    baselineParams.delta = max(double(fit_params.bump.delta), 0);
                else
                    baselineParams.delta = 0;
                end

                if isfield(fit_params, "optimo") && isstruct(fit_params.optimo)
                    baselineParams.optimo = fit_params.optimo;
                end

                % Best baseline plot method (keep actual A5 selection if present)
                bestBaselinePlot = bestBaseline;
                if ~ismember(bestBaselinePlot, methods)
                    bestBaselinePlot = local_choose_best_dynamic_baseline(validOverall, "bump_asymmetric", "bump_symmetric");
                end

                % Methods to include in the 5-subplot grid (fixed order)
                plotMethodsFixed = ["simple_selected", ...
                                    "optimo_lite", ...
                                    "bump_asymmetric", ...
                                    "optimo_lite_outcome_only", ...
                                    "bump_symmetric"];

                % Keep only those that exist in A5 outputs (but always keep simple_selected)
                plotMethods = strings(0,1);
                for i = 1:numel(plotMethodsFixed)
                    m = plotMethodsFixed(i);
                    if m == "simple_selected"
                        plotMethods(end+1) = m; %#ok<AGROW>
                    else
                        if ismember(m, methods)
                            plotMethods(end+1) = m; %#ok<AGROW>
                        end
                    end
                end

                % Build a per-participant wRMSE lookup from A5 tables (VALID split)
                wrmseLookup = local_build_participant_wrmse_lookup(A5plot, methods);

                % Iterate ALL validation participants (no selection)
                pid_list = string(read_participant_ids(participants_valid));
                allRows = table();
                allRows.participant_id = pid_list(:);
                writetable(allRows, fullfile(trajDir, "all_participants_for_trajectories.csv"));

                for k = 1:numel(participants_valid)
                    P = participants_valid(k);
                    pid = pid_list(k);

                    % Create participant subfolder
                    pidDir = fullfile(trajDir, char(pid));
                    if exist("ensure_dir","file") == 2
                        ensure_dir(pidDir);
                    else
                        if ~isfolder(pidDir), mkdir(pidDir); end
                    end

                    % Simulate model once to get aligned measurements & door events for overlays
                    simModel = trust_simulate_or_predict_one_participant("simple", theta_star, P, dt);
                    [t_meas, y_meas, kind_meas] = local_unpack_measurements(simModel.measurements);
                    doorEvents = [];
                    if isfield(simModel, "doorEvents"), doorEvents = simModel.doorEvents; end

                    % (A) Multi-panel grid (5 methods + note panel)
                    local_plot_methods_grid(pidDir, pid, plotMethods, ...
                        P, dt, weights, cfg, theta_star, baselineParams, ...
                        t_meas, y_meas, kind_meas, doorEvents, wrmseLookup);

                    % (B) Combo plot (constants dashed + model + best baseline)
                    constMethods = ["const_dispositional","const_global_train_mean","const_oracle_participant_mean"];
                    local_plot_combo(pidDir, pid, bestBaselinePlot, constMethods, ...
                        P, dt, weights, cfg, theta_star, baselineParams, ...
                        t_meas, y_meas, kind_meas, doorEvents, wrmseLookup);
                end
            end
        end
    end

    % ------------------------------------------------------------------
    % Terminal summary
    % ------------------------------------------------------------------
    fprintf("\n[Step A6] VALID overall wRMSE (lower is better):\n");
    fprintf("  %-30s | %-12s\n", "method", "wRMSE");
    fprintf("  %s\n", repmat('-',1,48));
    for i = 1:height(validOverall)
        fprintf("  %-30s | %-12.6g\n", string(validOverall.method(i)), validOverall.wRMSE_overall(i));
    end
    fprintf("\n[Step A6] Best baseline on VALID: %s (wRMSE=%.6g)\n", bestBaseline, bestVal);
    fprintf("[Step A6] Output: %s\n", outDir);
    fprintf("[Step A6] Trajectories: %s\n", trajDir);
end

% ---------------------------------------------------------------------
% Helper: choose best among two dynamic baselines using validOverall table
% ---------------------------------------------------------------------
function best = local_choose_best_dynamic_baseline(validOverall, m1, m2)
    best = m1;
    try
        if istable(validOverall) && ismember("method", validOverall.Properties.VariableNames) && ismember("wRMSE_overall", validOverall.Properties.VariableNames)
            T = validOverall;
            T.method = string(T.method);
            r1 = T.wRMSE_overall(T.method == m1);
            r2 = T.wRMSE_overall(T.method == m2);
            if ~isempty(r1) && ~isempty(r2) && isfinite(r2(1)) && isfinite(r1(1))
                if r2(1) < r1(1), best = m2; end
            elseif ~isempty(r2) && isfinite(r2(1))
                best = m2;
            end
        end
    catch
        best = m1;
    end
end

% ---------------------------------------------------------------------
% Build per-participant wRMSE lookup for VALID split from A5 tables
% Returns a containers.Map with key "method||participant_id" -> wRMSE double.
% ---------------------------------------------------------------------
function mp = local_build_participant_wrmse_lookup(A5plot, methods)
    mp = containers.Map('KeyType','char','ValueType','double');
    if ~isfield(A5plot, "allParticipants") || ~isstruct(A5plot.allParticipants)
        return;
    end

    for i = 1:numel(methods)
        m = string(methods(i));
        if ~(isfield(A5plot.allParticipants, m) && isfield(A5plot.allParticipants.(m), "valid"))
            continue;
        end
        Tp = A5plot.allParticipants.(m).valid;
        if ~istable(Tp) || ~all(ismember(["participant_id","wRMSE"], string(Tp.Properties.VariableNames)))
            continue;
        end

        pid = string(Tp.participant_id);
        wr  = double(Tp.wRMSE);

        for k = 1:numel(pid)
            key = char(m + "||" + pid(k));
            if isfinite(wr(k))
                mp(key) = wr(k);
            end
        end
    end
end

% ---------------------------------------------------------------------
% Helper: get per-participant wRMSE for a given method (NaN if missing)
% ---------------------------------------------------------------------
function v = local_lookup_participant_wrmse(wrmseLookup, method, pid)
    v = NaN;
    try
        key = char(string(method) + "||" + string(pid));
        if isKey(wrmseLookup, key)
            v = wrmseLookup(key);
        end
    catch
        v = NaN;
    end
end

% ---------------------------------------------------------------------
% NEW: Multi-panel grid plot (3x2) with 5 methods + note panel
%   - Panels 1..5: method trajectories (in given order)
%   - Panel 6 (bottom-right): axis off + wRMSE text box listing the 5 methods
% ---------------------------------------------------------------------
function local_plot_methods_grid(outDirPid, pid, plotMethods, P, dt, weights, cfg, theta_star, baselineParams, t_meas, y_meas, kind_meas, doorEvents, wrmseLookup)
    fig = figure('Visible','off','Color','w','Name',sprintf('Traj %s methods_grid', pid));
    local_make_fullscreen(fig);

    % Enforce up to 5 methods for panels 1..5
    plotMethods = string(plotMethods(:));
    if numel(plotMethods) > 5
        plotMethods = plotMethods(1:5);
    end

    % Plot each method in a dedicated subplot
    for i = 1:numel(plotMethods)
        method = plotMethods(i);
        subplot(3,2,i);
        hold on;

        sim = local_sim_for_plot(method, P, dt, weights, cfg, theta_star, baselineParams);
        plot(double(sim.t_grid(:)), double(sim.tau_hist(:)), '-', 'LineWidth', 1.8);

        local_plot_door_markers(doorEvents);
        local_plot_measurements_overlay(t_meas, y_meas, kind_meas);

        ylim([0 1]);
        grid on;
        xlabel('Time [s]');
        ylabel('Trust');

        wr = local_lookup_participant_wrmse(wrmseLookup, method, pid);
        if isfinite(wr)
            title(sprintf('%s  (wRMSE=%.4g)', char(method), wr), 'Interpreter','none');
        else
            title(sprintf('%s  (wRMSE=NA)', char(method)), 'Interpreter','none');
        end

        hold off;
    end

    % Bottom-right panel: wRMSE note box for quick scanning
    subplot(3,2,6);
    axis off;

    lines = strings(0,1);
    lines(end+1) = "Participant: " + string(pid); %#ok<AGROW>
    lines(end+1) = ""; %#ok<AGROW>
    lines(end+1) = "wRMSE (VALID):"; %#ok<AGROW>

    for i = 1:numel(plotMethods)
        method = plotMethods(i);
        wr = local_lookup_participant_wrmse(wrmseLookup, method, pid);
        lines(end+1) = sprintf("%d) %s : %s", i, char(method), char(local_fmt_num_or_na(wr))); %#ok<AGROW>
    end

    txt = strjoin(cellstr(lines), newline);
    local_add_note_in_empty_panel(txt);

    sgtitle(sprintf('Participant %s - methods (grid)', char(pid)), 'Interpreter','none');

    fname = sprintf("traj_%s_methods_grid.png", char(pid));
    exportgraphics(fig, fullfile(outDirPid, fname), 'Resolution', 220);
    close(fig);
end

% ---------------------------------------------------------------------
% Helper: plot combo (constants dashed + model + best baseline) + wRMSE notes
% ---------------------------------------------------------------------
function local_plot_combo(outDirPid, pid, bestBaseline, constMethods, P, dt, weights, cfg, theta_star, baselineParams, t_meas, y_meas, kind_meas, doorEvents, wrmseLookup)
    fig = figure('Visible','off','Color','w','Name',sprintf('Traj %s model_vs_best', pid));
    local_make_fullscreen(fig);
    hold on;

    % Constants (dashed)
    for i = 1:numel(constMethods)
        m = constMethods(i);
        sim = local_sim_for_plot(m, P, dt, weights, cfg, theta_star, baselineParams);
        plot(double(sim.t_grid(:)), double(sim.tau_hist(:)), '--', 'LineWidth', 1.2, 'DisplayName', m);
    end

    % Model (solid, thicker)
    simM = local_sim_for_plot("simple_selected", P, dt, weights, cfg, theta_star, baselineParams);
    plot(double(simM.t_grid(:)), double(simM.tau_hist(:)), '-', 'LineWidth', 2.2, 'DisplayName', "simple_selected");

    % Best baseline (solid)
    simB = local_sim_for_plot(bestBaseline, P, dt, weights, cfg, theta_star, baselineParams);
    plot(double(simB.t_grid(:)), double(simB.tau_hist(:)), '-', 'LineWidth', 2.0, 'DisplayName', bestBaseline);

    local_plot_door_markers(doorEvents);
    local_plot_measurements_overlay(t_meas, y_meas, kind_meas);

    ylim([0 1]);
    grid on;
    xlabel('Time [s]');
    ylabel('Trust');
    title(sprintf('Participant %s - model + best baseline + constants', pid), 'Interpreter','none');
    legend('Location','bestoutside');

    % ---- wRMSE annotation block for all plotted methods ----
    lines = strings(0,1);
    wrM = local_lookup_participant_wrmse(wrmseLookup, "simple_selected", pid);
    lines(end+1) = "simple_selected wRMSE: " + local_fmt_num_or_na(wrM); %#ok<AGROW>

    wrB = local_lookup_participant_wrmse(wrmseLookup, bestBaseline, pid);
    lines(end+1) = string(bestBaseline) + " wRMSE: " + local_fmt_num_or_na(wrB); %#ok<AGROW>

    for i = 1:numel(constMethods)
        m = constMethods(i);
        wrC = local_lookup_participant_wrmse(wrmseLookup, m, pid);
        lines(end+1) = string(m) + " wRMSE: " + local_fmt_num_or_na(wrC); %#ok<AGROW>
    end

    local_add_corner_note(strjoin(cellstr(lines), newline));

    hold off;

    fname = sprintf("traj_%s_model_vs_bestbaseline.png", char(pid));
    exportgraphics(fig, fullfile(outDirPid, fname), 'Resolution', 220);
    close(fig);
end

% ---------------------------------------------------------------------
% Helper: plot door event time markers (vertical dotted lines)
% ---------------------------------------------------------------------
function local_plot_door_markers(doorEvents)
    if isempty(doorEvents), return; end
    if ~isstruct(doorEvents) || ~isfield(doorEvents, "t"), return; end

    dts = [doorEvents.t];
    if isempty(dts), return; end
    yl = [0 1];
    for td = dts(:)'
        hh = plot([td td], yl, ':', 'LineWidth', 0.5);
        set(hh, 'HandleVisibility','off');
    end
end

% ---------------------------------------------------------------------
% Helper: run correct simulator for plotting
% ---------------------------------------------------------------------
function sim = local_sim_for_plot(method, P, dt, weights, cfg, theta_star, baselineParams)
    method = string(method);
    if method == "simple_selected"
        sim = trust_simulate_or_predict_one_participant("simple", theta_star, P, dt);
    else
        sim = trust_simulate_baseline_one_participant(method, P, dt, weights, baselineParams, cfg);
    end

    if ~isfield(sim, "t_grid") || ~isfield(sim, "tau_hist")
        error("[A6] Simulator output missing t_grid/tau_hist for method %s.", method);
    end
end

% ---------------------------------------------------------------------
% Helper: unpack aligned measurement struct array into vectors
% ---------------------------------------------------------------------
function [t_meas, y_meas, kind_meas] = local_unpack_measurements(measurements)
    if isempty(measurements)
        t_meas = zeros(0,1);
        y_meas = zeros(0,1);
        kind_meas = strings(0,1);
        return;
    end

    M = numel(measurements);
    t_meas    = NaN(M,1);
    y_meas    = NaN(M,1);
    kind_meas = strings(M,1);

    for i = 1:M
        if isfield(measurements(i), "t"),    t_meas(i) = double(measurements(i).t); end
        if isfield(measurements(i), "y"),    y_meas(i) = double(measurements(i).y); end
        if isfield(measurements(i), "kind"), kind_meas(i) = string(measurements(i).kind); end
    end

    ok = isfinite(t_meas) & isfinite(y_meas);
    t_meas = t_meas(ok);
    y_meas = y_meas(ok);
    kind_meas = kind_meas(ok);
end

% ---------------------------------------------------------------------
% Helper: rename table variable in a MATLAB-version-compatible way
% ---------------------------------------------------------------------
function T = local_rename_var(T, oldName, newName)
    oldName = char(oldName);
    newName = char(newName);
    vn = T.Properties.VariableNames;
    idx = find(strcmp(vn, oldName), 1);
    if isempty(idx)
        error("local_rename_var: variable '%s' not found.", oldName);
    end
    vn{idx} = newName;
    T.Properties.VariableNames = vn;
end

% ---------------------------------------------------------------------
% Helper: robust fullscreen sizing for trajectory figures
% ---------------------------------------------------------------------
function local_make_fullscreen(fig)
    try
        set(fig, 'Units','normalized');
        set(fig, 'Position',[0 0 1 1]);
    catch
    end
end

% ---------------------------------------------------------------------
% Overlay aligned measurements, grouped by kind (marker-only overlay)
% ---------------------------------------------------------------------
function local_plot_measurements_overlay(t_meas, y_meas, kind_meas)
    if isempty(t_meas) || isempty(y_meas)
        return;
    end

    t_meas = double(t_meas(:));
    y_meas = double(y_meas(:));
    kind_meas = string(kind_meas(:));

    kinds = unique(kind_meas);
    markers = {'o','s','^','d','x','+','v','>','<'};

    for i = 1:numel(kinds)
        k = kinds(i);
        mask = (kind_meas == k);
        mk = markers{1 + mod(i-1, numel(markers))};
        plot(t_meas(mask), y_meas(mask), mk, 'MarkerSize', 6, 'LineWidth', 1.0, 'DisplayName', "meas:" + k);
    end
end

% ---------------------------------------------------------------------
% Add a small text box in the bottom-right corner of the current axes
% ---------------------------------------------------------------------
function local_add_corner_note(txt)
    try
        ax = gca;
        x0 = ax.XLim(1);
        x1 = ax.XLim(2);
        y0 = ax.YLim(1);
        y1 = ax.YLim(2);

        % Place near bottom-right inside axes
        x = x1 - 0.02*(x1-x0);
        y = y0 + 0.02*(y1-y0);

        text(x, y, txt, ...
            'Units','data', ...
            'VerticalAlignment','bottom', ...
            'HorizontalAlignment','right', ...
            'FontSize', 10, ...
            'BackgroundColor','w', ...
            'EdgeColor',[0.3 0.3 0.3], ...
            'Margin', 6);
    catch
        % best-effort
    end
end

% ---------------------------------------------------------------------
% NEW: Add a text box inside an "empty" subplot panel (axes off)
% Uses normalized coordinates so it doesn't depend on axis limits.
% ---------------------------------------------------------------------
function local_add_note_in_empty_panel(txt)
    try
        axis off;
        text(0.98, 0.02, txt, ...
            'Units','normalized', ...
            'VerticalAlignment','bottom', ...
            'HorizontalAlignment','right', ...
            'FontSize', 11, ...
            'BackgroundColor','w', ...
            'EdgeColor',[0.3 0.3 0.3], ...
            'Margin', 8, ...
            'Interpreter','none');
    catch
        % best-effort
    end
end

function s = local_fmt_num_or_na(x)
    if ~isfinite(x)
        s = "NA";
    else
        s = sprintf("%.4g", x);
        s = string(s);
    end
end
