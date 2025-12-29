function params = trust_theta_to_params(theta)
% trust_theta_to_params  Map parameter vector to structured model parameters.
%
%   params = trust_theta_to_params(theta)
%
% This helper converts the global parameter vector theta (as used by
% fmincon / trust_cost_all) into a structured parameter representation
% used by the trust model components. No transformation or bounding of
% theta is performed here; this function is a pure reshaping/mapping
% utility.
%
% Parameter layout (current design; tau_offset is fixed, not optimized):
%   theta(1) = lambda_rep   (lambda_rep, reputation decay rate)        >= 0
%   theta(2) = phi_fail     (phi,  first failure magnitude)            in (0,1)
%   theta(3) = psi_succ     (psi,  first success magnitude)            in (0,1)
%   theta(4) = a_succ       (a,    success-shape parameter)            < 0
%   theta(5) = lambda_sit   (lambda_sit, situational risk sensitivity) > 0
%   theta(6) = lambda10     (lambda_10, base latent "above" rate)      > 0
%   theta(7) = kappa01      (kappa_01,  base latent "below" rate)      > 0
%
% Fixed design constants (not estimated from data):
%   params.lat.eps_lat       - deadzone around tau_disp for no drift
%   params.lat.gamma_above   - shaping parameter for "above" episodes
%   params.lat.epsilon_above - epsilon for "above" rate scaling
%   params.lat.gamma_below   - shaping parameter for "below" episodes
%   params.lat.epsilon_below - epsilon for "below" rate scaling
%   params.lat.tau_offset    - logistic offset for "below" dynamics
%
% These parameters are consumed by:
%   - trust_update_personal_experience (params.exp.*)
%   - trust_prepare_latent_sequence / trust_update_latent_sequence (params.lat.*)
%   - trust_update_reputation (params.rep.*)
%   - trust_compute_situational (params.sit.*)

    % Require at least the expected number of free parameters.
    if numel(theta) < 7
        error("trust_theta_to_params: expected at least 7 parameters, got %d.", numel(theta));
    end

    params = struct();

    % ------------------------------------------------------------
    % Personal experience (event-based updates at door trials)
    % ------------------------------------------------------------
    params.exp = struct();
    params.exp.phi_fail = theta(2);   % first failure magnitude
    params.exp.psi_succ = theta(3);   % first success magnitude
    params.exp.a_succ   = theta(4);   % success streak shape parameter

    % ------------------------------------------------------------
    % Latent trust dynamics (between door events)
    % ------------------------------------------------------------
    params.lat = struct();

    % Episode-specific base rates (estimated from data)
    params.lat.lambda10 = theta(6);   % base rate for "above" latent episodes
    params.lat.kappa01  = theta(7);   % base rate for "below" latent episodes

    % Small deadzone around tau_disp where no latent drift is applied.
    % This is kept fixed rather than estimated.
    params.lat.eps_lat = 1e-2;

    % Fixed hyperparameters for the lambda / kappa transformations.
    % These determine how the base rates are transformed to
    % episode-specific rates as a function of tau_disp and tau_lat0.
    % Defined by graphical inspection.
    params.lat.gamma_above   = 0.25;
    params.lat.epsilon_above = 0.02;

    params.lat.gamma_below   = 0.10;
    params.lat.epsilon_below = 0.01;
    params.lat.tau_offset = 0.1;

    % ------------------------------------------------------------
    % Reputation dynamics (global decay of initial reputation)
    % ------------------------------------------------------------
    params.rep = struct();
    params.rep.lambda_rep = theta(1);

    % ------------------------------------------------------------
    % Situational trust (instantaneous risk-based component)
    % ------------------------------------------------------------
    params.sit = struct();
    params.sit.lambda_sit = theta(5);

    % ------------------------------------------------------------
    % Dispositional trust
    % ------------------------------------------------------------
    % Dispositional trust is derived directly from questionnaires and
    % currently has no free parameters.
    params.disp = struct();
end
