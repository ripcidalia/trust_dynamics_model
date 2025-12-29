function stepM2_apply_14mapping(cleanMatPath, calibPath, outPath)
% stepM2_apply_14mapping  Map mid-block 14-item scores onto the 40-item scale.
%
%   stepM2_apply_14mapping(cleanMatPath, calibPath, outPath)
%
% This Step M2 script applies the global linear mapping obtained in Step M1
% to the mid-block 14-item questionnaire scores, converting them to their
% 40-item–equivalent percentages. The mapped values are added to each
% participant's questionnaires struct.
%
% In the broader pipeline:
%   - Step M1 (stepM1_14item_lopo) computes:
%        T40 ≈ a14 * T14 + b14
%     using LOPO cross-validation and global OLS.
%   - Step M2 uses the global coefficients (a14, b14) to map:
%        t14_mid1.total_percent        → t14_mid1.total_percent_40
%        t14_mid2.total_percent        → t14_mid2.total_percent_40
%     where all values are in percentage units (0..100).
%
% Inputs:
%   cleanMatPath : (optional) Path to MAT file containing 'participants_clean'
%                  (output of Step T1). Default:
%                  "derived/participants_time_stepT1.mat".
%
%   calibPath    : (optional) Path to MAT file containing 'calib14' from
%                  Step M1, with at least the fields:
%                      calib14.a14
%                      calib14.b14
%                  Default: "derived/measurement_step1_14item.mat".
%
%   outPath      : (optional) Output MAT file path.
%                  Default: "derived/participants_mapped14_stepM2.mat".
%
% Output (file):
%   Writes MAT file specified by outPath containing:
%
%     participants_mapped : struct array; a copy of participants_clean
%                           with additional fields:
%                               P(i).questionnaires.t14_mid1.total_percent_40
%                               P(i).questionnaires.t14_mid2.total_percent_40
%
%     info                : metadata struct with fields:
%                               .source_clean_file
%                               .calib_file
%                               .a14
%                               .b14
%                               .created
%                               .n_participants
%
% Assumptions:
%   - For each participant i, P(i).questionnaires.t14_mid1.total_percent
%     and t14_mid2.total_percent exist (or are missing/empty if the
%     respective questionnaire was not completed).
%   - These 14-item total_percent values are in 0..100.
%   - The mapping (a14, b14) was calibrated on compatible scales.

    if nargin < 1 || isempty(cleanMatPath)
        cleanMatPath = "derived/participants_time_stepT1.mat";
    end
    if nargin < 2 || isempty(calibPath)
        calibPath = "derived/measurement_step1_14item.mat";
    end
    if nargin < 3 || isempty(outPath)
        outPath = "derived/participants_mapped14_stepM2.mat";
    end

    if ~isfile(cleanMatPath)
        error("Clean participants file not found: %s", cleanMatPath);
    end
    if ~isfile(calibPath)
        error("Calibration file not found: %s", calibPath);
    end

    % ------------------------------------------------------------
    % 1) Load participants_clean (from Step 4)
    % ------------------------------------------------------------
    S = load(cleanMatPath, "participants_clean");
    if ~isfield(S, "participants_clean")
        error("File %s does not contain 'participants_clean'.", cleanMatPath);
    end
    participants = S.participants_clean;
    N = numel(participants);

    % ------------------------------------------------------------
    % 2) Load 14→40 calibration (a14, b14) from Step M1
    % ------------------------------------------------------------
    C = load(calibPath, "calib14");
    if ~isfield(C, "calib14")
        error("File %s does not contain 'calib14'.", calibPath);
    end
    calib14 = C.calib14;

    a14 = calib14.a14;
    b14 = calib14.b14;

    fprintf('[Step M2] Using 14→40 mapping: T40 = %.4f * T14 + %.4f\n', a14, b14);

    % ------------------------------------------------------------
    % 3) Apply linear mapping to t14_mid1 and t14_mid2 per participant
    % ------------------------------------------------------------
    for i = 1:N
        Q = participants(i).questionnaires;

        % --- t14_mid1 ---
        if isfield(Q, "t14_mid1") && ~isempty(Q.t14_mid1) ...
                && isfield(Q.t14_mid1, "total_percent")

            val = Q.t14_mid1.total_percent;
            val_num = coerce_scalar_double(val);
            if ~isnan(val_num)
                mapped = a14 * val_num + b14;
            else
                mapped = NaN;
            end
            participants(i).questionnaires.t14_mid1.total_percent_40 = mapped;
        else
            % Ensure field exists for consistency, even if NaN
            participants(i).questionnaires.t14_mid1.total_percent_40 = NaN;
        end

        % --- t14_mid2 ---
        if isfield(Q, "t14_mid2") && ~isempty(Q.t14_mid2) ...
                && isfield(Q.t14_mid2, "total_percent")

            val = Q.t14_mid2.total_percent;
            val_num = coerce_scalar_double(val);
            if ~isnan(val_num)
                mapped = a14 * val_num + b14;
            else
                mapped = NaN;
            end
            participants(i).questionnaires.t14_mid2.total_percent_40 = mapped;
        else
            participants(i).questionnaires.t14_mid2.total_percent_40 = NaN;
        end
    end

    % ------------------------------------------------------------
    % 4) Save updated participants with mapped mid-block scores
    % ------------------------------------------------------------
    if ~isfolder(fileparts(outPath))
        mkdir(fileparts(outPath));
    end

    participants_mapped = participants;
    info = struct();
    info.source_clean_file = cleanMatPath;
    info.calib_file        = calibPath;
    info.a14               = a14;
    info.b14               = b14;
    info.created           = char(datetime('now','Format','yyyy-MM-dd HH:mm:ss'));
    info.n_participants    = N;

    save(outPath, "participants_mapped", "info", "-v7.3");

    fprintf('[Step M2] Applied 14→40 mapping to mid-block questionnaires for %d participants.\n', N);
    fprintf('          Saved to %s\n', outPath);
end

% -------------------------------------------------------------------------
% Local helper: coerce_scalar_double
% -------------------------------------------------------------------------
function v = coerce_scalar_double(val)
% COERCE_SCALAR_DOUBLE  Convert numeric/string/cell to scalar double or NaN.
%
%   v = coerce_scalar_double(val)
%
% Attempts to coerce "val" into a scalar double:
%   - numeric: returns double(val) if scalar, otherwise NaN
%   - string/char: strips non-numeric characters and uses str2double
%   - cell: unwraps the first element and retries
% Returns NaN if conversion is not possible.

    v = NaN;

    if iscell(val)
        if isempty(val), return; end
        val = val{1};
    end

    if isnumeric(val)
        if isscalar(val)
            v = double(val);
        end
        return;
    end

    if isstring(val) || ischar(val)
        s = regexprep(strtrim(string(val)), '[^0-9\.\-eE]+', '');
        vv = str2double(s);
        if isscalar(vv)
            v = vv;
        end
        return;
    end
end
