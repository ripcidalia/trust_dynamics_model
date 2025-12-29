function summary = analyze_participants_overview()
% analyze_participants_overview  Participant-level overview plots and tables.
%
% This script loads the time-enriched participants struct array and produces:
%   - Demographics distributions (age range, gender)
%   - Set (door-sequence) counts (SetA–SetH)
%   - Review expected distribution
%   - Expected vs response review scatter (all participants)
%   - Questionnaire score distributions (overall + by set)
%   - Device type and browser summary tables (for reporting)
%   - Emergency choice summary table (for reporting)
%   - Balance matrices (columns = sets; rows = category bins) for:
%       age range, gender, 40pre bins, review expected bins, emergency choice
%
% Output:
%   summary  - struct containing extracted vectors and key summary tables.
%
% Files saved:
%   figs/participant_analysis/*.png and *.fig
%   tables/participant_analysis/*.csv
%   derived/participant_analysis.mat

    %% ------------------------------------------------------------
    % 0) Paths and load
    % -------------------------------------------------------------
    inMat = "derived/participants_time_stepT1.mat";
    if ~isfile(inMat)
        error("File not found: %s. Run preprocessing (Step T1) first.", inMat);
    end

    S = load(inMat, "participants_clean");
    if ~isfield(S, "participants_clean")
        error("Variable participants_clean not found in %s.", inMat);
    end
    participants = S.participants_clean;
    N = numel(participants);

    % Use existing top-level folders (create if missing)
    figsRoot = "figs";
    tabsRoot = "tables";
    if ~isfolder(figsRoot), mkdir(figsRoot); end
    if ~isfolder(tabsRoot), mkdir(tabsRoot); end

    % Keep outputs organized without nesting under derived/
    outFigs = fullfile(figsRoot, "participant_analysis");
    outTabs = fullfile(tabsRoot, "participant_analysis");
    if ~isfolder(outFigs), mkdir(outFigs); end
    if ~isfolder(outTabs), mkdir(outTabs); end

    outMat = "derived/participant_analysis.mat";

    fprintf('[analyze_participants_overview] Loaded %d participants from %s\n', N, inMat);

    %% ------------------------------------------------------------
    % 1) Extract fields (assumed complete/consistent after preprocessing)
    % -------------------------------------------------------------
    pid      = strings(N,1);
    set_id   = strings(N,1);
    device   = strings(N,1);
    browser  = strings(N,1);

    age_rng  = strings(N,1);
    gender   = strings(N,1);

    rev_exp  = nan(N,1);
    rev_resp = nan(N,1);

    emerg    = strings(N,1);

    q40pre   = nan(N,1);
    q40post  = nan(N,1);
    q14m1    = nan(N,1);
    q14m2    = nan(N,1);

    for i = 1:N
        P = participants(i);

        pid(i)    = string(P.participant_id);
        set_id(i) = string(P.set_id);
        device(i) = string(P.device_type);
        browser(i)= string(P.browser_name);

        age_rng(i)= string(P.demographics.age_range);
        gender(i) = string(P.demographics.gender);

        % Reviews: nested struct/cell access as provided.
        rev_exp(i)  = double(P.reviews.items(1).response_struct{1,1}.expected);
        rev_resp(i) = double(P.reviews.items(1).response_struct{1,1}.response);

        % Emergency: "self" or "robot"
        emerg(i) = string(P.emergency.choice);

        % Questionnaires: numeric scalar or struct with total_percent
        q40pre(i)  = local_get_questionnaire_total_percent(P.questionnaires.t40_pre);
        q40post(i) = local_get_questionnaire_total_percent(P.questionnaires.t40_post);
        q14m1(i)   = local_get_questionnaire_total_percent(P.questionnaires.t14_mid1);
        q14m2(i)   = local_get_questionnaire_total_percent(P.questionnaires.t14_mid2);
    end

    % Normalize set IDs to a consistent ordering (SetA..SetH where possible)
    set_id = local_normalize_set_id(set_id);

    % Fixed set order for consistent plotting (SetA..SetH). All participants
    % are assumed to belong to these sets after preprocessing.
    setCats = "Set" + string(('A':'H')');
    if ~all(ismember(set_id, setCats))
        bad = unique(set_id(~ismember(set_id, setCats)));
        error(['Unexpected set_id values found (expected only SetA..SetH). ' ...
               'Found: %s'], strjoin(bad, ", "));
    end
    setC = categorical(set_id, setCats, 'Ordinal', true);

    %% ------------------------------------------------------------
    % 2) Basic tables for reporting
    % -------------------------------------------------------------
    T_device = local_count_table(device, "device_type");
    writetable(T_device, fullfile(outTabs, "participants_by_device_type.csv"));

    T_browser = local_count_table(browser, "browser_name");
    writetable(T_browser, fullfile(outTabs, "participants_by_browser.csv"));

    T_emerg = local_count_table(emerg, "emergency_choice");
    writetable(T_emerg, fullfile(outTabs, "participants_by_emergency_choice.csv"));

    %% ------------------------------------------------------------
    % 3) Plots
    % -------------------------------------------------------------
    local_bar_counts(age_rng, "Age range", ...
        "Participants by age range", fullfile(outFigs, "age_range"));

    local_bar_counts(gender, "Gender", ...
        "Participants by gender", fullfile(outFigs, "gender"));

    % Participants by set
    f = figure('Name','Participants by set','Color','w');
    countsBySet = countcats(setC);
    bar(countsBySet);
    grid on;
    set(gca,'XTick',1:numel(setCats),'XTickLabel',cellstr(setCats));
    xlabel('Door-trial sequence set');
    ylabel('Number of participants');
    title('Participants by door-trial sequence set');
    local_save_figure(f, fullfile(outFigs, "participants_by_set"));

    % Expected review distribution
    f = figure('Name','Expected review distribution','Color','w');
    histogram(rev_exp);
    grid on;
    xlabel('Expected review score');
    ylabel('Count');
    title('Distribution of expected review scores');
    local_save_figure(f, fullfile(outFigs, "review_expected_hist"));

    % Expected vs response review (scatter + y=x)
    f = figure('Name','Expected vs response review','Color','w');
    scatter(rev_exp, rev_resp, 30, 'filled');
    grid on;
    xlabel('Expected review score');
    ylabel('Response review score');
    title('Expected vs response review (all participants)');

    xMin = min([rev_exp; rev_resp]);
    xMax = max([rev_exp; rev_resp]);
    hold on;
    plot([xMin xMax], [xMin xMax], 'k--', 'LineWidth', 1);
    hold off;

    local_save_figure(f, fullfile(outFigs, "review_expected_vs_response_scatter"));

    % Questionnaires (overall histogram + by-set boxplot)
    local_questionnaire_plots(q40pre,  setC, '40pre',  outFigs);
    local_questionnaire_plots(q40post, setC, '40post', outFigs);
    local_questionnaire_plots(q14m1,   setC, '14mid1', outFigs);
    local_questionnaire_plots(q14m2,   setC, '14mid2', outFigs);

    %% ------------------------------------------------------------
    % 4) Balance matrices (columns = sets; rows = category bins)
    % -------------------------------------------------------------
    % 4.1 Age range x set
    local_balance_heatmap( ...
        categorical(age_rng), setC, ...
        "age_range", "set_id", ...
        "Balance matrix: age range vs set", ...
        fullfile(outFigs, "balance_age_range_vs_set"), ...
        fullfile(outTabs, "balance_age_range_vs_set.csv"));

    % 4.2 Gender x set (robust to any strings, including self-describe)
    local_balance_heatmap( ...
        categorical(gender), setC, ...
        "gender", "set_id", ...
        "Balance matrix: gender vs set", ...
        fullfile(outFigs, "balance_gender_vs_set"), ...
        fullfile(outTabs, "balance_gender_vs_set.csv"));

    % 4.3 40pre bins x set (fixed bins, questionnaires are 0–100)
    [q40pre_bin, q40pre_bin_labels] = local_bin_percent(q40pre);
    local_balance_heatmap( ...
        categorical(q40pre_bin, q40pre_bin_labels, 'Ordinal', true), setC, ...
        "q40pre_bin", "set_id", ...
        "Balance matrix: 40pre (binned) vs set", ...
        fullfile(outFigs, "balance_40pre_bins_vs_set"), ...
        fullfile(outTabs, "balance_40pre_bins_vs_set.csv"));

    % 4.4 Review expected (robust: discrete numeric values if few uniques,
    % otherwise falls back to quantile bins)
    [revexp_bin, revexp_bin_labels] = local_bin_numeric_discrete_or_quantiles(rev_exp);
    local_balance_heatmap( ...
        categorical(revexp_bin, revexp_bin_labels, 'Ordinal', true), setC, ...
        "review_expected_bin", "set_id", ...
        "Balance matrix: expected review vs set", ...
        fullfile(outFigs, "balance_review_expected_vs_set"), ...
        fullfile(outTabs, "balance_review_expected_vs_set.csv"));

    % 4.5 Emergency choice x set
    local_balance_heatmap( ...
        categorical(emerg), setC, ...
        "emergency_choice", "set_id", ...
        "Balance matrix: emergency choice vs set", ...
        fullfile(outFigs, "balance_emergency_vs_set"), ...
        fullfile(outTabs, "balance_emergency_vs_set.csv"));

    %% ------------------------------------------------------------
    % 5) Save consolidated MAT in derived/
    % -------------------------------------------------------------
    summary = struct();
    summary.N = N;

    summary.participant_id = pid;
    summary.set_id = set_id;
    summary.device_type = device;
    summary.browser_name = browser;

    summary.age_range = age_rng;
    summary.gender = gender;
    summary.emergency_choice = emerg;

    summary.review_expected = rev_exp;
    summary.review_response = rev_resp;

    summary.q40pre_percent  = q40pre;
    summary.q40post_percent = q40post;
    summary.q14mid1_percent = q14m1;
    summary.q14mid2_percent = q14m2;

    summary.tables.device = T_device;
    summary.tables.browser = T_browser;
    summary.tables.emergency = T_emerg;

    save(outMat, "summary");

    fprintf('[analyze_participants_overview] Done.\n');
    fprintf('  Figures: %s\n', outFigs);
    fprintf('  Tables : %s\n', outTabs);
    fprintf('  MAT    : %s\n', outMat);

end

%% =====================================================================
%  Helper functions (local)
% ======================================================================

function v = local_get_questionnaire_total_percent(Q)
% Supports either:
%   - numeric scalar (already total percent)
%   - struct with field total_percent
%   - struct with nested fields (fallback tries a few common names)
    if isnumeric(Q) && isscalar(Q)
        v = double(Q);
        return;
    end
    if isstruct(Q)
        if isfield(Q, 'total_percent')
            v = double(Q.total_percent);
            return;
        end
        if isfield(Q, 'total') && isnumeric(Q.total) && isscalar(Q.total)
            v = double(Q.total);
            return;
        end
        if isfield(Q, 'score') && isnumeric(Q.score) && isscalar(Q.score)
            v = double(Q.score);
            return;
        end
    end
    error('Questionnaire value could not be read as a scalar percent.');
end

function set_id = local_normalize_set_id(set_id)
% Ensures values like "A", "SetA", "set a" become "SetA" where possible.
    set_id = string(set_id);
    set_id = strip(set_id);

    for i = 1:numel(set_id)
        s = strip(set_id(i));
        s = replace(s, " ", "");
        sUp = upper(s);

        % Strip leading "SET" if present
        if startsWith(sUp, "SET")
            sUp = extractAfter(sUp, 3);
        end

        % If single letter A-H, normalize to SetX
        if strlength(sUp) == 1 && sUp >= "A" && sUp <= "H"
            set_id(i) = "Set" + sUp;
            continue;
        end

        % If ends with A-H, also normalize (e.g., "SetA", "Condition_SetB")
        lastChar = extractAfter(sUp, strlength(sUp)-1);
        if strlength(lastChar) == 1 && lastChar >= "A" && lastChar <= "H"
            set_id(i) = "Set" + lastChar;
        else
            % Unknown/other: keep original (trimmed)
            set_id(i) = s;
        end
    end
end

function T = local_count_table(strVec, varName)
% Count occurrences of a string/categorical vector into a table.
    c = categorical(string(strVec));
    cats = categories(c);
    counts = countcats(c);

    % VariableNames must be a cell array of character vectors (R2022b).
    vn1 = char(string(varName));

    T = table(string(cats), counts, ...
        'VariableNames', {vn1, 'count'});
    T = sortrows(T, 'count', 'descend');
end

function local_bar_counts(strVec, xLabelStr, titleStr, outBase)
% Bar chart for categorical counts.
    c = categorical(string(strVec));
    cats = categories(c);
    counts = countcats(c);

    f = figure('Name', titleStr, 'Color', 'w');
    bar(counts);
    grid on;
    set(gca,'XTick',1:numel(cats),'XTickLabel',cellstr(cats));
    xlabel(xLabelStr);
    ylabel('Number of participants');
    title(titleStr);

    local_save_figure(f, outBase);
end

function local_questionnaire_plots(qPercent, setC, tag, outFigs)
% Produce:
%   - overall histogram of questionnaire total_percent
%   - boxplot grouped by set
    f = figure('Name', sprintf('%s overall histogram', tag), 'Color', 'w');
    histogram(qPercent);
    grid on;
    xlabel(sprintf('%s total (percent)', tag));
    ylabel('Count');
    title(sprintf('Distribution of %s questionnaire scores (all participants)', tag));
    local_save_figure(f, fullfile(outFigs, sprintf('%s_hist_all', tag)));

    f = figure('Name', sprintf('%s by set', tag), 'Color', 'w');
    boxplot(qPercent, setC);
    grid on;
    xlabel('Door-trial sequence set');
    ylabel(sprintf('%s total (percent)', tag));
    title(sprintf('%s questionnaire scores by set', tag));
    local_save_figure(f, fullfile(outFigs, sprintf('%s_box_by_set', tag)));
end

function local_balance_heatmap(rowCat, colCat, rowName, colName, titleStr, outFigBase, outCsvPath)
% Build counts matrix and show as heatmap; also export as CSV.
    rowCat = categorical(rowCat);
    colCat = categorical(colCat);

    rowLevels = categories(rowCat);
    colLevels = categories(colCat);

    M = zeros(numel(rowLevels), numel(colLevels));
    for i = 1:numel(rowLevels)
        for j = 1:numel(colLevels)
            M(i,j) = sum(rowCat == rowLevels{i} & colCat == colLevels{j});
        end
    end

    % Export table
    T = array2table(M, 'VariableNames', matlab.lang.makeValidName(colLevels));
    T = addvars(T, string(rowLevels), 'Before', 1, 'NewVariableNames', rowName);
    writetable(T, outCsvPath);

    % Heatmap figure
    f = figure('Name', titleStr, 'Color', 'w');
    h = heatmap(colLevels, rowLevels, M);
    h.XLabel = colName;
    h.YLabel = rowName;
    h.Title  = titleStr;

    local_save_figure(f, outFigBase);
end

function [binStr, binLabels] = local_bin_percent(x)
% Bin a percent-valued vector into fixed ranges.
% Bins: [0,20), [20,40), [40,60), [60,80), [80,100]
    edges = [0 20 40 60 80 100.0001];
    labels = ["0–20", "20–40", "40–60", "60–80", "80–100"];

    binIdx = discretize(x, edges);
    binStr = strings(size(x));
    for k = 1:numel(x)
        binStr(k) = labels(binIdx(k));
    end
    binLabels = labels;
end

function [binStr, binLabels] = local_bin_numeric_discrete_or_quantiles(x)
% If x looks discrete (few unique finite values), bin by unique values.
% Otherwise bin into 5 quantile-based bins.
    x = x(:);
    xFinite = x(isfinite(x));
    u = unique(xFinite);

    % "Discrete" heuristic: small number of unique values and close to integers
    isNearlyInteger = all(abs(u - round(u)) < 1e-9);
    if numel(u) <= 12 && isNearlyInteger
        u = sort(u);
        binLabels = string(u);
        binStr = strings(size(x));
        for k = 1:numel(x)
            binStr(k) = string(x(k));
        end
        return;
    end

    % Quantile bins (5)
    q = quantile(xFinite, [0 0.2 0.4 0.6 0.8 1.0]);
    q(1) = q(1) - 1e-9;
    q(end) = q(end) + 1e-9;

    edges = unique(q);
    if numel(edges) < 3
        edges = linspace(min(xFinite), max(xFinite) + 1e-9, 6);
    end

    idx = discretize(x, edges);
    nBins = max(idx);

    binLabels = strings(nBins,1);
    for i = 1:nBins
        a = edges(i);
        b = edges(i+1);
        binLabels(i) = sprintf('%.3g–%.3g', a, b);
    end

    binStr = strings(size(x));
    for k = 1:numel(x)
        binStr(k) = binLabels(idx(k));
    end
end

function local_save_figure(f, outBase)
% Save .png and .fig with a consistent base path.
    pngPath = outBase + ".png";
    figPath = outBase + ".fig";
    saveas(f, pngPath);
    savefig(f, figPath);
    close(f);
end
