% main.m
%
% MAIN  Entry point for the SAR trust + behavior modelling pipeline.
%
% This script runs the end-to-end workflow for the human–robot trust study,
% from raw data processing through trust model fitting, behavioral modelling,
% robustness analyses, and coupled rollout validation.
%
% -------------------------------------------------------------------------
% Pipeline overview
%
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
%   C) Analysis run creation + reporting
%      (run-local outputs under derived/analysis_runs/<run_id>/)
%
%      4)  stepA0_import_fit_results(cfg)
%          - Register optimisation outputs and create an analysis run ID
%
%      5)  stepA1_prepare_analysis_run(run_id)
%          - Archive inputs for this run (participant splits, mappings, metadata)
%
%      6)  stepA2_simple_mode_fit_report(run_id, resultsMatPath)
%          - Generate trust model fit quality summaries and figures
%
%      7)  stepA3_select_theta_for_coupled(run_id, resultsMatPath)
%          - Select theta_star for downstream prediction and coupled analyses
%
%      8)  stepA4_sensitivity_simple_mode(run_id, resultsMatPath)
%          - Sensitivity analyses around fitted trust parameters
%
%      9)  stepA5_compare_baselines_simple_mode(run_id, resultsMatPath)
%          - Compare trust model baselines using held-out evaluation metrics
%
%     10)  stepA6_report_baseline_comparison_simple_mode(run_id)
%          - Compile baseline-comparison reporting artifacts
%
%     11)  stepA7_build_behavior_dataset(run_id)
%          - Build per-door behavioral dataset (follow/override decisions)
%
%     12)  stepA8_behavior_fit_eval(run_id)
%          - Fit behavioral models on TRAIN and evaluate on VALID
%
%     13)  stepA9_behavior_rollouts(run_id)
%          - Coupled generative rollout analysis using global behavior models
%
%     14)  stepA10_behavior_fit_by_participant(run_id)
%          - Fit behavioral models per participant (door-resampling bootstrap)
%          - Select best model via BIC
%
%     15)  stepA11_behavior_param_robustness(run_id)
%          - Robustness analysis of behavioral parameters
%            * random door splits (quantitative)
%            * blockwise stability (qualitative: stable / drifting / under-identified)
%
%     16)  stepA12_behavior_rollouts_personalized(run_id)
%          - Personalized coupled rollouts using per-participant behavior params
%          - Guardrail fallback for under-identified or invalid fits
%
%     17)  stepA13_trust_divergence_sanity_check(run_id)
%          - Trust realism diagnostics:
%            * divergence between simple replay and coupled trust trajectories
%            * effect of personalization vs global behavior parameters
%            * full-time-grid divergence metrics and visualizations
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
stepA7_build_behavior_dataset(run_id, "Split", "train");
stepA7_build_behavior_dataset(run_id, "Split", "valid");
stepA8_behavior_fit_eval(run_id);
stepA9_behavior_rollouts(run_id);
stepA10_behavior_fit_by_participant(run_id);
stepA11_behavior_param_robustness(run_id);
stepA12_behavior_rollouts_personalized(run_id);
stepA13_trust_divergence_sanity_check(run_id);
