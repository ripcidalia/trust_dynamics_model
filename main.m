% ExperimentalDataProcessing/main.m
%
% MAIN  Entry point for the SAR trust modelling pipeline.
%
% This script runs the full data-processing pipeline for the human–robot
% trust study:
%   1) Preprocessing of raw experimental events into cleaned, time-stamped
%      participant data and derived intermediate files.
%   2) Postprocessing and calibration of trust measurements, including
%      questionnaire and probe mappings and measurement weights.
%
% The actual logic is implemented in:
%   - preprocessing.m : Steps 1–4 and T1 (data loading, cleaning, timing).
%   - postprocessing.m: Steps M1–M6 (measurement mappings and diagnostics).
%
% This script is intended to be run from the project root.

clear; clc;
addpath(genpath("src"));
preprocessing
postprocessing


method = "ga";
dt = 1;
preset = "overnight";
theta0 = [ % Initial conditions for parameters
    3e-3
    0.15
    0.10
    0.20
    1.0
    1e-4
    1e-4
];

[theta_hat_ga, fval_ga, exitflag_ga, output_ga] = fit_trust_parameters(method, dt, preset, theta0);

save("results.mat", "theta_hat_ga", "fval_ga", "exitflag_ga", "output_ga");

% cfg = struct();
% cfg.dt = 1.0; % Simulation time step
% cfg.theta0 = [ % Initial conditions for parameters
%     3e-3
%     0.15
%     0.10
%     0.20
%     1.0
%     1e-4
%     1e-4
% ];
% cfg.run_tag = "train1";
% cfg.results_dir = "derived/fit_runs";
% cfg.checkpoint_dir = "derived/checkpoints";
% 
% results = run_trust_optimisation_pipeline(cfg);

