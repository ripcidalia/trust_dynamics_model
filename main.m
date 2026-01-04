% main.m
%
% MAIN  Entry point for the SAR trust + behavior modelling pipeline.
%
% This script runs the end-to-end workflow for the human–robot trust study,
% from raw data processing through parameter fitting and analysis reporting.
%
% Pipeline overview
%   A) Data processing (project-wide, produces cleaned participant structs)
%      1) preprocessing.m
%         - Load raw experimental logs/events
%         - Clean/standardize event streams and timestamps
%         - Produce per-participant, time-aligned data suitable for modelling
%      2) postprocessing.m
%         - Map questionnaire/probe responses into model-ready measurements
%         - Apply measurement weights/calibration and generate diagnostics
%
%   B) Trust dynamics model fitting (produces fit artifacts in derived/fit_runs)
%      3) run_trust_optimisation_pipeline(cfg)
%         - Fit trust-dynamics parameters (theta) in "simple" mode using dt
%         - Save run results and checkpoints for reproducibility
%
%   C) Analysis run creation + reporting (run-local, under derived/analysis_runs/<run_id>/)
%      4) stepA0_import_fit_results(cfg)
%         - Register the optimisation outputs and create an analysis run id
%      5) stepA1_prepare_analysis_run(run_id)
%         - Archive inputs for this run (participants splits, mappings, metadata)
%      6) stepA2_simple_mode_fit_report(run_id, resultsMatPath)
%         - Generate fit quality summaries/figures for the selected run
%      7) stepA3_select_theta_for_coupled(run_id, resultsMatPath)
%         - Select theta_star used for downstream prediction/coupled analyses
%      8) stepA4_sensitivity_simple_mode(run_id, resultsMatPath)
%         - Sensitivity analyses around the fitted parameters
%      9) stepA5_compare_baselines_simple_mode(run_id, resultsMatPath)
%         - Compare trust model baselines and produce evaluation metrics
%     10) stepA6_report_baseline_comparison_simple_mode(run_id)
%         - Compile baseline-comparison reporting artifacts
%
% Notes
%   - Steps A7–A9 (behavior dataset, behavioral model fit/eval, coupled rollouts)
%     are not invoked in this script; they are typically run after A1/A3
%     using the created run_id.
%   - This script is intended to be executed from the project root so that
%     relative paths under "derived/" resolve correctly.
%
% -------------------------------------------------------------------------

clear; clc;
addpath(genpath("src"));
preprocessing
postprocessing

cfg = struct();
cfg.dt = 1.0;  % Simulation time step
cfg.theta0 = [ % Initial conditions for parameters
    3e-3       % lambda_rep
    1          % alpha_sit
    1          % lambda_sit
    0.5        % phi_fail
    0.45       % phi_succ
    0.10       % a_succ
    1e-4       % lambda_lat
    1e-4       % kappa_lat
];
cfg.run_tag = "thesis_fit";
cfg.results_dir = "derived/fit_runs";
cfg.checkpoint_dir = "derived/checkpoints";

results = run_trust_optimisation_pipeline(cfg);

[run_id, resultsMatPath] = stepA0_import_fit_results(cfg);
stepA1_prepare_analysis_run(run_id);
stepA2_simple_mode_fit_report(run_id, resultsMatPath);
stepA3_select_theta_for_coupled(run_id, resultsMatPath);
stepA4_sensitivity_simple_mode(run_id, resultsMatPath);
stepA5_compare_baselines_simple_mode(run_id, resultsMatPath);
stepA6_report_baseline_comparison_simple_mode(run_id);
stepA7_build_behavior_dataset(run_id);
stepA8_behavior_fit_eval(run_id);
stepA9_behavior_rollouts(run_id);
