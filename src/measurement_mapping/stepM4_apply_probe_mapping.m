function stepM4_apply_probe_mapping(cleanMatPath, calibProbePath, outPath)
% stepM4_apply_probe_mapping  Map all trust probes onto the 40-item scale.
%
%   stepM4_apply_probe_mapping(cleanMatPath, calibProbePath, outPath)
%
% This Step M4 function applies the global probe→40 mapping obtained in
% Step M3 to every single trust probe for all participants. It appends
% a new field, value_40, to each probe:
%
%       value_40 = a1 * value + b1
%
% where "value" is the original probe score (0..100) and (a1, b1) are the
% global probe→40 coefficients estimated in stepM3_probe_calibration.
%
% In the broader pipeline:
%   - Step M3 produces calibProbe.a1 and calibProbe.b1 (global map),
%     and an associated probe variance and weight used later in WLS.
%   - Step M4 uses (a1, b1) to precompute 40-scale equivalents for all
%     probes, so that later stages (e.g., trust_cost_* functions) can
%     consume trust measurements on a unified scale.
%
% Inputs:
%   cleanMatPath   : (optional) Path to MAT file containing participants_mapped,
%                    typically the output of Step M2. Default:
%                        "derived/participants_mapped14_stepM2.mat"
%
%   calibProbePath : (optional) Path to MAT file containing calibProbe from
%                    Step M3, with fields:
%                        calibProbe.a1
%                        calibProbe.b1
%                    Default:
%                        "derived/measurement_step3_probe.mat"
%
%   outPath        : (optional) Output MAT file path.
%                    Default: "derived/participants_mapped14_stepM2.mat".
%
% Output (file):
%   Writes MAT file specified by outPath containing:
%       participants_probes_mapped : participants_mapped with each element
%           P(i).trustProbes(k).value_40 populated using the global map.
%       info                       : struct with metadata (source paths,
%                                    mapping coefficients, timestamps).
%
% Assumptions:
%   - participants_mapped(i).trustProbes is an array of probe structs.
%   - Each probe struct has a field "value" containing a numeric score
%     (0..100) or NaN.
%   - No structure or field names are modified; only value_40 is added.

    if nargin < 1 || isempty(cleanMatPath)
        cleanMatPath = "derived/participants_mapped14_stepM2.mat";
    end
    if nargin < 2 || isempty(calibProbePath)
        calibProbePath = "derived/measurement_step3_probe.mat";
    end
    if nargin < 3 || isempty(outPath)
        outPath = "derived/participants_probes_mapped_stepM4.mat";
    end

    if ~isfile(cleanMatPath)
        error("participants_mapped file not found: %s", cleanMatPath);
    end
    if ~isfile(calibProbePath)
        error("calibProbe file not found: %s", calibProbePath);
    end

    % ------------------------------------------------------------
    % 1) Load participants and probe calibration
    % ------------------------------------------------------------
    S = load(cleanMatPath, "participants_mapped");
    participants = S.participants_mapped;

    C = load(calibProbePath, "calibProbe");
    calibProbe = C.calibProbe;
    a1 = calibProbe.a1;
    b1 = calibProbe.b1;

    % ------------------------------------------------------------
    % 2) Apply probe→40 mapping to all probes for all participants
    % ------------------------------------------------------------
    N = numel(participants);
    for i = 1:N
        U = participants(i).trustProbes;
        for k = 1:numel(U)
            val = U(k).value;
            type = string(U(k).questionnaire_type);
            if ~isnan(val)
                if type == "t40_pre"
                    participants(i).trustProbes(k).value_40 = participants(i).questionnaires.t40_pre.total_percent;
                    participants(i).trustProbes(k).t_s = participants(i).questionnaires.t40_pre.t_s;                    
                elseif type == "t40_post"
                    participants(i).trustProbes(k).value_40 = participants(i).questionnaires.t40_post.total_percent;
                    participants(i).trustProbes(k).t_s = participants(i).questionnaires.t40_post.t_s;
                elseif type == "t14_mid1"
                    participants(i).trustProbes(k).value_40 = participants(i).questionnaires.t14_mid1.total_percent_40;
                    participants(i).trustProbes(k).t_s = participants(i).questionnaires.t14_mid1.t_s;
                elseif type == "t14_mid2"
                    participants(i).trustProbes(k).value_40 = participants(i).questionnaires.t14_mid2.total_percent_40;
                    participants(i).trustProbes(k).t_s = participants(i).questionnaires.t14_mid2.t_s;
                else
                    % Map raw probe value (0..100) to 40-item equivalent
                    participants(i).trustProbes(k).value_40 = a1 * double(val) + b1;
                end
            else
                % Preserve NaN for missing or invalid probe values
                participants(i).trustProbes(k).value_40 = NaN;
            end        
        end
    end

    % ------------------------------------------------------------
    % 3) Save updated participants with mapped probes
    % ------------------------------------------------------------
    if ~isfolder(fileparts(outPath))
        mkdir(fileparts(outPath));
    end

    participants_probes_mapped = participants;

    info = struct();
    info.source_mapped_file = cleanMatPath;
    info.calib_probe_file  = calibProbePath;
    info.a1                = a1;
    info.b1                = b1;
    info.created           = char(datetime('now','Format','yyyy-MM-dd HH:mm:ss'));
    info.n_participants    = N;

    save(outPath, "participants_probes_mapped", "info", "-v7.3");

    fprintf('[Step M4] Applied probe→40 mapping to all probes for %d participants.\n', N);
    fprintf('          Saved to %s\n', outPath);
end
