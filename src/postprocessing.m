% ExperimentalDataProcessing/postprocessing.m
%
% POSTPROCESSING  Run the calibration and measurement-mapping pipeline.
%
% This script orchestrates the postprocessing steps (M1–M6) for the
% human–robot trust model in the SAR study. It operates on the
% time-annotated participant data produced by the preprocessing pipeline
% and generates:
%   - Mappings between the 14-item and 40-item trust questionnaires.
%   - Calibrated mappings from single-probe trust measurements to the
%     40-item scale.
%   - Measurement weights for different trust instruments.
%   - Reputation-related diagnostics and summary figures.
%
% The steps are:
%   M1: Leave-one-participant-out (LOPO) mapping of 14-item to 40-item scores.
%   M2: Apply the 14-item→40-item mapping to mid-block 14-item scores.
%   M3: Calibrate probe→40 mapping (Phase A + Phase B).
%   M4: Apply the probe→40 mapping to all probes.
%   M5: Compute and save measurement weights.
%   M6: Run reputation bias diagnostics and generate associated plots/PDFs.
%
% This script is intended to be run from the project root and assumes:
%   - Time-annotated participants at derived/participants_time_stepT1.mat.
%   - Derived measurement outputs are written under derived/.
%   - Figures and reputation diagnostics are written under figs/ (M6).
%
% No computational logic is implemented here; this file only sequences the
% calls to dedicated functions in src/. The function calls and file paths
% should remain aligned with the rest of the pipeline.

%% Step M1: 14-item LOPO cross-validation and global mapping
clear; clc; addpath(genpath("src"));

% Runs LOPO cross-validation to derive a mapping from the mid-block
% 14-item questionnaire scores to the 40-item scale and saves the results.
stepM1_14item_lopo("derived/participants_train_stepV.mat");

%% Step M2: Apply 14-item→40-item mapping to mid-block scores
clear; clc; addpath(genpath("src"));

% Applies the mapping derived in M1 to mid-block 14-item scores, producing
% 40-item-equivalent trust levels for those time points.
stepM2_apply_14mapping("derived/participants_train_stepV.mat", ...
                       "derived/measurement_step1_14item.mat");

%% Step M3: Probe calibration (Phase A + Phase B)

% Calibrates the mapping from single trust probes to the 40-item scale
% using the 14-item→40-item mapping as reference.
stepM3_probe_calibration("derived/participants_train_stepV.mat", ...
                         "derived/measurement_step1_14item.mat");

%% Step M4: Apply probe→40 mapping to all probes

% Uses the calibrated probe→40 mapping (from M3) to convert all trust
% probe measurements to the common 40-item-equivalent scale.
stepM4_apply_probe_mapping("derived/participants_mapped14_stepM2.mat", ...
                           "derived/measurement_step3_probe.mat");

%% Step M5: Compute and save measurement weights

% Computes measurement weights (e.g., w40, w14, w_probe) that will later
% be used in the trust cost function and parameter estimation.
stepM5_save_measurement_weights( ...
    "derived/measurement_step1_14item.mat", ...
    "derived/measurement_step3_probe.mat");

%% Step M6: Reputation diagnostics and figure generation

% Evaluates reputation-related effects and bias, and generates diagnostic
% plots and PDFs stored under the specified directory.
stepM6_reputation_bias("derived/participants_train_stepV.mat", "figs");
