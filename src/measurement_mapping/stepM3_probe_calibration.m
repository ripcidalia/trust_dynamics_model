function stepM3_probe_calibration(cleanMatPath, calib14Path)
% stepM3_probe_calibration  Calibrate single-probe trust measurements.
%
%   stepM3_probe_calibration(cleanMatPath, calib14Path)
%
% This Step M3 function performs a two-phase calibration of the single
% trust probes (slider ratings) onto the 40-item questionnaire scale.
% It uses both pre/post questionnaire anchors and mid-block 14-item
% questionnaires to estimate the measurement variance of the probes.
%
% In the broader pipeline:
%   - Step M1 (14-item LOPO) provides:
%         calib14.sigma2_14   : variance of 14-item mapping error
%         calib14.a14_lopo    : participant-specific 14→40 slopes
%         calib14.b14_lopo    : participant-specific 14→40 intercepts
%   - Step M3 estimates:
%         a1, b1              : global probe→40 linear map
%         sigma2_1A           : probe variance from pre/post anchors
%         sigma2_1B           : additional probe variance from midpoints
%         sigma2_probe        : combined probe variance
%         w_probe             : probe weight for WLS = 1 / sigma2_probe
%
% Phase A (anchors):
%   Uses only pre/post 40-item questionnaires and the probes taken
%   immediately after these questionnaires to:
%       1) estimate a global linear mapping probe → 40-item,
%       2) estimate a Phase-A probe variance sigma2_1A via LOPO.
%
% Phase B (midpoints):
%   Uses mid-block 14-item questionnaires and their associated probes
%   (after_questionnaire origin) together with participant-specific LOPO
%   14→40 mappings to obtain midpoint residuals and decompose:
%       Var( midpoint residual ) = sigma2_14 + sigma2_1B
%   which gives an estimate of additional probe variance sigma2_1B.
%
% The final probe variance is a weighted combination of Phase A and Phase B:
%       sigma2_probe = (nA * sigma2_1A + nB * sigma2_1B) / (nA + nB)
% and the probe weight is:
%       w_probe = 1 / sigma2_probe
%
% Inputs:
%   cleanMatPath : (optional) Path to MAT file containing participants_clean
%                  (output of Step T1). Default:
%                  "derived/participants_time_stepT1.mat".
%
%   calib14Path  : (optional) Path to MAT file containing calib14 from
%                  stepM1_14item_lopo, with at least fields:
%                      calib14.sigma2_14
%                      calib14.a14_lopo
%                      calib14.b14_lopo
%                  Default:
%                  "derived/measurement_step1_14item.mat".
%
% Output (file):
%   Writes "derived/measurement_step3_probe.mat" containing struct
%   calibProbe with fields:
%       .a1, .b1              global probe→40 mapping coefficients
%       .sigma2_1A            Phase-A probe variance from anchors
%       .sigma2_1B            additional Phase-B probe variance
%       .sigma2_probe         combined probe variance
%       .w_probe              probe measurement weight (1 / sigma2_probe)
%       .resA_pre             pre-anchor LOPO residuals
%       .resA_post            post-anchor LOPO residuals
%       .ids                  participant IDs
%       .R_B                  midpoint residuals (14 vs probe)
%       .nA, .nB              counts of residuals used in A and B
%
% Assumptions:
%   - For each participant:
%       * P.questionnaires.t40_pre.total_percent and
%         P.questionnaires.t40_post.total_percent are available (0..100).
%       * P.trustProbes is an array where:
%             .origin == "after_questionnaire"
%             .questionnaire_type is "t40_pre", "t40_post",
%                                   "t14_mid1", or "t14_mid2".
%             .value is the probe reading (0..100).
%       * P.questionnaires.t14_mid1.total_percent and
%         P.questionnaires.t14_mid2.total_percent exist when the mid-block
%         questionnaires are present.
%   - calib14 was computed on the same participant set and aligns with
%     participants_clean in order.

    if nargin < 1 || isempty(cleanMatPath)
        cleanMatPath = "derived/participants_time_stepT1.mat";
    end
    if nargin < 2 || isempty(calib14Path)
        calib14Path = "derived/measurement_step1_14item.mat";
    end

    if ~isfile(cleanMatPath)
        error("participants_clean file not found: %s", cleanMatPath);
    end
    if ~isfile(calib14Path)
        error("calib14 file not found: %s", calib14Path);
    end

    % ------------------------------------------------------------
    % 1) Load participants and 14-item calibration
    % ------------------------------------------------------------
    % Load participants_clean from preprocessing
    S = load(cleanMatPath, "participants_clean");
    participants = S.participants_clean;
    N = numel(participants);

    % Load calib14 (14-item LOPO results)
    C = load(calib14Path, "calib14");
    calib14   = C.calib14;
    sigma2_14 = calib14.sigma2_14;
    a14_lopo  = calib14.a14_lopo;
    b14_lopo  = calib14.b14_lopo;

    % ------------------------------------------------------------
    % 2) Phase A: anchors (pre/post questionnaires and their probes)
    % ------------------------------------------------------------
    % For each participant we extract:
    %   Q40_pre, Q40_post    : questionnaire scores (0..100)
    %   P_pre,  P_post       : probe values after these questionnaires

    Q40_pre   = nan(N,1);
    Q40_post  = nan(N,1);
    P_pre     = nan(N,1);
    P_post    = nan(N,1);
    ids       = strings(N,1);

    for i = 1:N
        Pi     = participants(i);
        ids(i) = string(Pi.participant_id);

        q = Pi.questionnaires;
        Q40_pre(i)  = double(q.t40_pre.total_percent);
        Q40_post(i) = double(q.t40_post.total_percent);

        % Locate probes whose origin is after the corresponding questionnaire
        P_pre(i)  = find_probe_after_q(Pi.trustProbes, "t40_pre");
        P_post(i) = find_probe_after_q(Pi.trustProbes, "t40_post");
    end

    % LOPO on anchors to estimate probe variance sigma2_1A
    resA_pre  = nan(N,1);
    resA_post = nan(N,1);

    for i = 1:N
        % Training set: all other participants
        J = setdiff(1:N, i);

        x_train = [P_pre(J);  P_post(J)];    % probe values (anchors)
        y_train = [Q40_pre(J); Q40_post(J)]; % 40-item targets

        % Simple linear regression: y = a1_i * x + b1_i
        A     = [x_train, ones(numel(x_train),1)];
        theta = A \ y_train;
        a1_i  = theta(1);
        b1_i  = theta(2);

        % Predict held-out participant's anchors and form residuals
        yhat_pre  = a1_i * P_pre(i)  + b1_i;
        yhat_post = a1_i * P_post(i) + b1_i;

        resA_pre(i)  = Q40_pre(i)  - yhat_pre;
        resA_post(i) = Q40_post(i) - yhat_post;
    end

    % Pooled LOPO variance over all anchor residuals
    E_A = [resA_pre; resA_post];
    E_A = E_A(~isnan(E_A));
    if numel(E_A) < 2
        error("Not enough non-NaN anchor residuals to estimate probe Phase-A variance.");
    end
    eA_bar    = mean(E_A);
    sigma2_1A = sum((E_A - eA_bar).^2) / (numel(E_A) - 1);
    nA        = numel(E_A);

    % Global probe→40 mapping (no LOPO), used in Phase B
    Xall        = [P_pre; P_post];
    Yall        = [Q40_pre; Q40_post];
    Aall        = [Xall, ones(numel(Xall),1)];
    theta_glob  = Aall \ Yall;
    a1          = theta_glob(1);
    b1          = theta_glob(2);

    % ------------------------------------------------------------
    % 3) Phase B: midpoints (variance decomposition with 14-item LOPO map)
    % ------------------------------------------------------------
    % For each midpoint m=1,2 of participant i:
    %   - Q14_mid(i,m) is mapped to 40-scale using per-participant LOPO
    %     (a14_lopo(i), b14_lopo(i)).
    %   - P_mid(i,m) is mapped to 40-scale using global probe→40 (a1, b1).
    %   - Residual r = Q40_14_hat - Q40_probe_hat is stored in R_B.

    R_B = [];  % accumulate midpoint residuals across participants

    for i = 1:N
        Pi = participants(i);

        % Mid-block 14-item scores (if present)
        q       = Pi.questionnaires;
        hasMid1 = isfield(q,"t14_mid1") && ~isempty(q.t14_mid1) && isfield(q.t14_mid1,"total_percent");
        hasMid2 = isfield(q,"t14_mid2") && ~isempty(q.t14_mid2) && isfield(q.t14_mid2,"total_percent");

        % Participant-specific LOPO 14→40 mapping
        a14_i = a14_lopo(i);
        b14_i = b14_lopo(i);
        if isnan(a14_i) || isnan(b14_i)
            % If LOPO mapping is not defined for this participant,
            % skip them in Phase B.
            continue;
        end

        % --- midpoint 1 ---
        if hasMid1
            Q14_mid1 = coerce_scalar_double(q.t14_mid1.total_percent);
            P_mid1   = find_probe_after_q(Pi.trustProbes, "t14_mid1");

            if ~isnan(Q14_mid1) && ~isnan(P_mid1)
                % Map 14-item midpoint to 40-scale using LOPO map
                Q40_14_hat = a14_i * Q14_mid1 + b14_i;
                % Map probe midpoint to 40-scale using global probe→40
                Q40_P_hat  = a1    * P_mid1   + b1;
                % Residual between these two 40-equivalent estimates
                R_B(end+1,1) = Q40_14_hat - Q40_P_hat; %#ok<AGROW>
            end
        end

        % --- midpoint 2 ---
        if hasMid2
            Q14_mid2 = coerce_scalar_double(q.t14_mid2.total_percent);
            P_mid2   = find_probe_after_q(Pi.trustProbes, "t14_mid2");

            if ~isnan(Q14_mid2) && ~isnan(P_mid2)
                Q40_14_hat = a14_i * Q14_mid2 + b14_i;
                Q40_P_hat  = a1    * P_mid2   + b1;
                R_B(end+1,1) = Q40_14_hat - Q40_P_hat; %#ok<AGROW>
            end
        end
    end

    % Estimate additional probe variance sigma2_1B from midpoint residuals
    if isempty(R_B)
        warning("No midpoint pairs found for probe Phase B; using Phase A variance only.");
        sigma2_1B = 0;
        nB        = 0;
    else
        R_B = R_B(~isnan(R_B));
        if numel(R_B) >= 2
            rB_bar       = mean(R_B);
            sigma2_delta = sum((R_B - rB_bar).^2) / (numel(R_B) - 1);
            % Variance decomposition:
            %   Var(R_B) = sigma2_14 + sigma2_1B
            sigma2_1B = max(sigma2_delta - sigma2_14, 0);
            nB        = numel(R_B);
        else
            sigma2_1B = 0;
            nB        = numel(R_B);
        end
    end

    % ------------------------------------------------------------
    % 4) Combine Phase A and B into final probe variance and weight
    % ------------------------------------------------------------
    if nA + nB > 0
        sigma2_probe = (nA * sigma2_1A + nB * sigma2_1B) / (nA + nB);
    else
        error("No residuals available to estimate overall probe variance.");
    end
    w_probe = 1 / sigma2_probe;

    % ------------------------------------------------------------
    % 5) Pack results and save to disk
    % ------------------------------------------------------------
    calibProbe = struct();
    calibProbe.a1           = a1;
    calibProbe.b1           = b1;
    calibProbe.sigma2_1A    = sigma2_1A;
    calibProbe.sigma2_1B    = sigma2_1B;
    calibProbe.sigma2_probe = sigma2_probe;
    calibProbe.w_probe      = w_probe;
    calibProbe.resA_pre     = resA_pre;
    calibProbe.resA_post    = resA_post;
    calibProbe.ids          = ids;
    calibProbe.R_B          = R_B;
    calibProbe.nA           = nA;
    calibProbe.nB           = nB;

    if ~isfolder("derived")
        mkdir("derived");
    end
    outPath = "derived/measurement_step3_probe.mat";
    save(outPath, "calibProbe", "-v7.3");

    fprintf('[Step M3] Probe calibration completed on %d participants.\n', N);
    fprintf('          sigma2_1A = %.4f, sigma2_1B = %.4f, sigma2_probe = %.4f, w_probe = %.4f\n', ...
            sigma2_1A, sigma2_1B, sigma2_probe, w_probe);
end

% -------------------------------------------------------------------------
% Local helpers
% -------------------------------------------------------------------------

function v = find_probe_after_q(trustProbes, qtype)
% find_probe_after_q  Extract probe value following a given questionnaire.
%
%   v = find_probe_after_q(trustProbes, qtype)
%
% Returns the numeric probe value (0..100) for the first probe whose
% origin is "after_questionnaire" and whose questionnaire_type equals
% the requested qtype (e.g., "t40_pre", "t40_post", "t14_mid1", "t14_mid2").
% Returns NaN if no such probe is found or if value is missing.

    v = NaN;
    for k = 1:numel(trustProbes)
        tp = trustProbes(k);
        if isfield(tp,"origin") && tp.origin == "after_questionnaire" ...
                && isfield(tp,"questionnaire_type") && tp.questionnaire_type == qtype
            if isfield(tp,"value")
                v = double(tp.value);
                return;
            end
        end
    end
end

function v = coerce_scalar_double(val)
% coerce_scalar_double  Convert numeric/string/cell to scalar double or NaN.
%
%   v = coerce_scalar_double(val)
%
% Used here for questionnaire percentage values that may be stored as
% numeric, string, or cell. Returns NaN if conversion fails.

    v = NaN;
    if iscell(val)
        if isempty(val), return; end
        val = val{1};
    end
    if isnumeric(val)
        if isscalar(val), v = double(val); end
        return;
    end
    if isstring(val) || ischar(val)
        s  = regexprep(strtrim(string(val)), '[^0-9\.\-eE]+', '');
        vv = str2double(s);
        if isscalar(vv), v = vv; end
        return;
    end
end
