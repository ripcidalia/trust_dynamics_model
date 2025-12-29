function stepM1_14item_lopo(cleanMatPath)
% stepM1_14item_lopo  LOPO weighting and global mapping for 14-item scores.
%
%   stepM1_14item_lopo(cleanMatPath)
%
% This Step M1 script:
%   1) Loads cleaned participant data from Step 4.
%   2) Extracts, for each participant, the pre/post 40-item questionnaire
%      scores and the corresponding 14-item–equivalent scores derived from
%      the 40-item questionnaire (anchors).
%   3) Performs leave-one-participant-out (LOPO) cross-validation to fit a
%      linear mapping:
%           Q40 ≈ a14 * Q14_eq + b14
%      where Q14_eq is the 14-item–equivalent score derived from the
%      40-item questionnaire.
%   4) Uses LOPO residuals to estimate a pooled residual variance
%      sigma2_14 and measurement weight w14 = 1 / sigma2_14.
%   5) Computes a global (non-LOPO) OLS 14→40 mapping over all anchors.
%   6) Saves all calibration results to:
%         derived/measurement_step1_14item.mat
%
% In the broader pipeline:
%   - The global mapping (a14, b14) is used later to transform mid-block
%     14-item questionnaire scores onto the 40-item scale.
%   - The pooled LOPO variance sigma2_14 defines the relative weight of
%     these measurements in the overall WLS cost.
%
% Inputs:
%   cleanMatPath : (optional) Path to a MAT file containing the variable
%                  'participants_clean' (output of Step T1).
%                  Default: "derived/participants_time_stepT1.mat".
%
% Outputs (file):
%   The MAT file derived/measurement_step1_14item.mat containing struct:
%
%   calib14 : struct with fields
%       .a14        - global OLS slope for Q14_eq → Q40 mapping
%       .b14        - global OLS intercept for Q14_eq → Q40 mapping
%       .sigma2_14  - pooled LOPO residual variance across all participants
%       .w14        - measurement weight 1 / sigma2_14
%       .res_pre    - N×1 vector of LOPO residuals at 40-pre anchors
%       .res_post   - N×1 vector of LOPO residuals at 40-post anchors
%       .ids        - N×1 string array of participant IDs
%       .a14_lopo   - N×1 LOPO slopes (fit with each participant left out)
%       .b14_lopo   - N×1 LOPO intercepts
%
% Assumptions:
%   - Each participant P(i) has a field P(i).questionnaires with subfields
%       t40_pre.total_percent
%       t40_post.total_percent
%       t40_pre.trust14_equiv_total_percent
%       t40_post.trust14_equiv_total_percent
%     all in units of percentage (0..100).
%   - There are at least 3 participants, so the LOPO regression is well
%     defined for each held-out participant.

    if nargin < 1 || isempty(cleanMatPath)
        cleanMatPath = "derived/participants_time_stepT1.mat";
    end
    if ~isfile(cleanMatPath)
        error("Clean participants file not found: %s", cleanMatPath);
    end

    % ------------------------------------------------------------
    % 1) Load cleaned participants from Step 4
    % ------------------------------------------------------------
    S = load(cleanMatPath, "participants_clean");
    if ~isfield(S, "participants_clean")
        error("File %s does not contain 'participants_clean'.", cleanMatPath);
    end
    P = S.participants_clean;
    N = numel(P);
    if N < 3
        error("Need at least 3 participants for meaningful LOPO; found %d.", N);
    end

    % ------------------------------------------------------------
    % 2) Extract anchor pairs (Q14_eq, Q40) at pre and post
    %
    % For each participant i we build:
    %   Q40_pre(i),  Q40_post(i)  : true 40-item percentages
    %   Q14_pre(i),  Q14_post(i)  : 14-item–equivalent percentages derived
    %                              from the 40-item questionnaire.
    % ------------------------------------------------------------
    Q40_pre  = nan(N,1);
    Q40_post = nan(N,1);
    Q14_pre  = nan(N,1);  % 14-equivalent at pre
    Q14_post = nan(N,1);  % 14-equivalent at post
    ids      = strings(N,1);

    for i = 1:N
        ids(i) = string(P(i).participant_id);

        q = P(i).questionnaires;

        Q40_pre(i)  = double(q.t40_pre.total_percent);
        Q40_post(i) = double(q.t40_post.total_percent);

        Q14_pre(i)  = double(q.t40_pre.trust14_equiv_total_percent);
        Q14_post(i) = double(q.t40_post.trust14_equiv_total_percent);
    end

    % ------------------------------------------------------------
    % 3) LOPO cross-validation to estimate sigma2_14
    %
    % For each held-out participant i:
    %   - Fit linear mapping on all others J = setdiff(1:N,i):
    %         y ≈ a14_i * x + b14_i
    %   - Predict Q40_pre(i) and Q40_post(i) from their 14-equivalent
    %     counterparts.
    %   - Store residuals for later pooling.
    % ------------------------------------------------------------
    res_pre   = nan(N,1);
    res_post  = nan(N,1);
    a14_lopo  = nan(N,1);
    b14_lopo  = nan(N,1);

    for i = 1:N
        J = setdiff(1:N, i);

        % Training data: stacked pre and post anchors
        x_train = [Q14_pre(J); Q14_post(J)];
        y_train = [Q40_pre(J); Q40_post(J)];

        % Linear regression: [a14_i; b14_i] = argmin ||A*theta - y||
        A = [x_train, ones(numel(x_train),1)];

        theta = A \ y_train;
        a14_i = theta(1);
        b14_i = theta(2);

        a14_lopo(i) = a14_i;
        b14_lopo(i) = b14_i;

        % Predict held-out anchors for participant i
        yhat_pre  = a14_i * Q14_pre(i)  + b14_i;
        yhat_post = a14_i * Q14_post(i) + b14_i;

        res_pre(i)  = Q40_pre(i)  - yhat_pre;
        res_post(i) = Q40_post(i) - yhat_post;
    end

    % ------------------------------------------------------------
    % 4) Pooled LOPO residual variance and weight
    %
    % We pool pre and post residuals to obtain sigma2_14 and define the
    % measurement weight w14 = 1 / sigma2_14.
    % ------------------------------------------------------------
    E = [res_pre; res_post];
    E = E(~isnan(E));
    if numel(E) < 2
        error("Not enough non-NaN residuals to estimate variance.");
    end
    e_bar     = mean(E);
    sigma2_14 = sum((E - e_bar).^2) / (numel(E) - 1);
    w14       = 1 / sigma2_14;

    % ------------------------------------------------------------
    % 5) Global 14→40 OLS mapping (no LOPO)
    %
    % This uses all anchor pairs simultaneously to estimate a single
    % linear map:
    %       Q40 ≈ a14 * Q14_eq + b14
    % ------------------------------------------------------------
    Xall = [Q14_pre; Q14_post];
    Yall = [Q40_pre; Q40_post];
    Aall = [Xall, ones(numel(Xall),1)];
    theta_glob = Aall \ Yall;

    a14 = theta_glob(1);
    b14 = theta_glob(2);

    % ------------------------------------------------------------
    % 6) Pack results and save to disk
    % ------------------------------------------------------------
    calib14 = struct();
    calib14.a14        = a14;
    calib14.b14        = b14;
    calib14.sigma2_14  = sigma2_14;
    calib14.w14        = w14;
    calib14.res_pre    = res_pre;
    calib14.res_post   = res_post;
    calib14.ids        = ids;
    calib14.a14_lopo   = a14_lopo;   % per-participant LOPO slope
    calib14.b14_lopo   = b14_lopo;   % per-participant LOPO intercept

    if ~isfolder("derived")
        mkdir("derived");
    end

    outPath = "derived/measurement_step1_14item.mat";
    save(outPath, "calib14", "-v7.3");

    fprintf('[Step M1] 14-item LOPO completed on %d participants.\n', N);
    fprintf('          sigma2_14 = %.4f, w14 = %.4f, a14 = %.4f, b14 = %.4f\n', ...
            sigma2_14, w14, a14, b14);
end
