function stepM5_save_measurement_weights(calib14Path, calibProbePath)
% stepM5_save_measurement_weights  Consolidate measurement weights for all instruments.
%
%   stepM5_save_measurement_weights(calib14Path, calibProbePath)
%
% Step M5 creates a unified measurement-weights structure that combines:
%   - 14-item questionnaire variance/weight (from Step M1),
%   - Probe variance/weight (from Step M3),
%   - A reference weight for the 40-item questionnaire.
%
% The resulting weights struct is used in the trust-fitting cost functions
% to build a Weighted Least Squares (WLS) objective that balances the three
% measurement types on a common variance scale.
%
% Inputs (optional):
%   calib14Path    - Path to MAT file containing 'calib14' from
%                    stepM1_14item_lopo. Default:
%                        "derived/measurement_step1_14item.mat"
%
%   calibProbePath - Path to MAT file containing 'calibProbe' from
%                    stepM3_probe_calibration. Default:
%                        "derived/measurement_step3_probe.mat"
%
% Output (file):
%   Writes "derived/measurement_weights.mat" containing:
%       weights.w40           - reference weight for 40-item scores (set to 1)
%       weights.w14           - weight for 14-item scores (1 / sigma2_14)
%       weights.w_probe       - weight for probes (1 / sigma2_probe)
%       weights.sigma2_14     - LOPO-estimated variance of 14-item mapping
%       weights.sigma2_probe  - combined variance of probe mapping
%       weights.info          - metadata (source files, timestamp, description)
%
% Assumptions:
%   - calib14Path MAT file contains a struct 'calib14' with fields:
%         calib14.w14, calib14.sigma2_14
%   - calibProbePath MAT file contains a struct 'calibProbe' with fields:
%         calibProbe.w_probe, calibProbe.sigma2_probe
%   - No fields are modified; this step only packages and saves them.

    % -----------------------------
    % Handle file inputs
    % -----------------------------
    if nargin < 1 || isempty(calib14Path)
        calib14Path = "derived/measurement_step1_14item.mat";
    end
    if nargin < 2 || isempty(calibProbePath)
        calibProbePath = "derived/measurement_step3_probe.mat";
    end

    if ~isfile(calib14Path)
        error("14-item calibration file not found: %s", calib14Path);
    end
    if ~isfile(calibProbePath)
        error("Probe calibration file not found: %s", calibProbePath);
    end

    % -----------------------------
    % Load calibration results
    % -----------------------------
    C14 = load(calib14Path, "calib14");
    CPr = load(calibProbePath, "calibProbe");

    calib14    = C14.calib14;
    calibProbe = CPr.calibProbe;

    % -----------------------------
    % Build merged weights struct
    % -----------------------------
    % 40-item questionnaire is treated as the reference measurement with
    % unit weight; other instruments are down/up-weighted relative to this.
    weights = struct();
    weights.w40          = 1;                     % reference scale
    weights.w14          = calib14.w14;           % 1 / sigma2_14 (14-item)
    weights.w_probe      = calibProbe.w_probe;    % 1 / sigma2_probe (probes)

    weights.sigma2_14    = calib14.sigma2_14;
    weights.sigma2_probe = calibProbe.sigma2_probe;

    % Metadata about how these weights were constructed
    weights.info = struct();
    weights.info.calib14_source    = calib14Path;
    weights.info.calibProbe_source = calibProbePath;
    weights.info.created           = char(datetime('now', ...
                                        'Format','yyyy-MM-dd HH:mm:ss'));
    weights.info.description = ...
        "Unified measurement weights for 40-item, 14-item, and probe instruments.";

    % -----------------------------
    % Save
    % -----------------------------
    if ~isfolder("derived")
        mkdir("derived");
    end

    outPath = "derived/measurement_weights.mat";
    save(outPath, "weights", "-v7.3");

    fprintf("[Step M5] Measurement weights saved to %s\n", outPath);
    fprintf("          w40 = %.3f, w14 = %.3f, w_probe = %.3f\n", ...
        weights.w40, weights.w14, weights.w_probe);
end
