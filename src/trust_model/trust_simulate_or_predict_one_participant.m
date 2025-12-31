function sim = trust_simulate_or_predict_one_participant(mode, theta, P, dt, steepness)
% trust_simulate_or_predict_one_participant  Simulate or predict trust dynamics for one participant.
%
%   sim = trust_simulate_or_predict_one_participant(mode, theta, P, dt, steepness)
%
% This function runs a full forward simulation of the trust model for a
% single participant, given a parameter vector θ and a time step size.
% It builds or uses a cached version of the time grid, attaches all door
% events and trust measurements, steps the trust dynamics forward in time,
% and records:
%   - the trajectory of total trust and its components, and
%   - model predictions at the times where trust is measured.
%
% The function will either use the recorded door event trials outcomes to deterministically
% update trust dynamics (mode: "simple"), or couple the behavioral model to predict user
% choices at door trials (mode: "coupled"). In the case of coupling with the behavioral model,
% another parameter ("steepness") can be passed to control the steepness of the probabilistic
% transition between the two behavioral states (follow/not follow).
%
% Inputs:
%   mode       Simulation mode (string):
%                - "simple", or "coupled"
%
%   theta      Parameter vector (column or row) that is interpreted by
%              trust_theta_to_params. The exact layout of theta is handled
%              entirely inside trust_theta_to_params.
%
%   P          Participant struct, typically coming from:
%                - participants_probes_mapped_stepM4.mat
%              and enriched by:
%                - stepT1_add_times (P.doorTrials(k).t_s, P.trustProbes(j).t_s, etc.)
%              It is assumed that build_time_grid_and_events(P, dt) succeeds
%              and that trust_init_state can initialize a valid state.
%
%   dt         Time step in seconds for the simulation time grid.
%              If omitted or empty, dt defaults to 1 second.
%
%   steepness  Controls the steepness of the probabilistic transition between the two
%              behavioral actions (follow/ not follow). Only relevant four "coupled" mode.
%              If omitted or empty, steepness defaults to 10.
%
% Output:
%   sim        Struct with simulation outputs:
%                 .t_grid        [K x 1] vector of time points (seconds)
%                 .doorEvents    [nDoor x 1] struct array (from
%                                get_cached_time_grid_and_events)
%                 .tau_hist      [K x 1] total trust trajectory τ(t)
%                 .tau_lat_hist  [K x 1] latent component trajectory τ_lat(t)
%                 .tau_rep_hist  [K x 1] reputation component trajectory τ_rep(t)
%                 .tau_sit_hist  [K x 1] situational component trajectory τ_sit(t)
%                 .risk_hist     [K x 1] risk value history r(t)
%                 .measurements  struct array (from build_time_grid_and_events),
%                                including measurement times and ground-truth y
%                 .y_hat         [nMeas x 1] model predictions at those
%                                measurement times, aligned with
%                                sim.measurements(m).y
%
% Notes:
%   - The evolution within each time step is performed by trust_step, which
%     updates all trust components and internal substates.
%   - Multiple door events are allowed at the same grid time t_k and are
%     applied sequentially.
%   - Measurements (questionnaires and probes) do not alter the trust state;
%     they only induce a read-out of the current τ at their respective times.

    if isempty(mode)
        warning("trust_simulate_or_predict_one_participant: no mode specified. Proceeding with simple mode (no behavioral prediction).");
        mode = "simple";
    end

    if mode == "coupled" && (nargin < 5 || isempty(steepness))
        warning("trust_simulate_or_predict_one_participant: no steepness specified for coupled mode. Proceeding with default steepness = 10.");
        steepness = 10;
    end
    
    if nargin < 4 || isempty(dt)
        dt = 1;
    end

    mode = string(lower(mode));

    % ------------------------------------------------------------
    % 1) Map θ → params
    % ------------------------------------------------------------
    params = trust_theta_to_params(theta);

    % ------------------------------------------------------------
    % 2) Build time grid and attach door / measurement events
    %    (CACHED: depends only on participant data + dt, not on theta)
    % ------------------------------------------------------------
    [t_grid, doors_at_k, doorEvents, meas_at_k, measurements] = ...
        get_cached_time_grid_and_events(P, dt);

    K     = numel(t_grid);
    nMeas = numel(measurements);

    % ------------------------------------------------------------
    % 3) Initialize trust state at t ≈ 0 (pre-40 anchor)
    % ------------------------------------------------------------
    state = trust_init_state(params, P);

    % Ensure latent component, last_risk and self-confidence are present.
    % These should normally be set by trust_init_state; the checks below
    % serve as a safety net in case older participants or data formats are
    % encountered.
    if ~isfield(state, "tau_lat") || isempty(state.tau_lat)
        % Fallback: latent = total - reputation (ignoring situational)
        state.tau_lat = state.tau - state.tau_rep;
    end
    if ~isfield(state, "last_risk")
        state.last_risk = NaN;
    end
    if ~isfield(state, "sc")
        % trust_init_state should set sc from emergency choice; this is a final fallback
        state.sc = 0.5;
    end

    % Latent-episode substate (used by trust_update_latent_sequence)
    state.lat = struct();
    state.lat.mode       = "none";
    state.lat.lambda_seq = NaN;
    state.lat.kappa_seq  = NaN;
    state.lat.sigma      = NaN;
    state.lat.tau0       = state.tau_lat;

    % ------------------------------------------------------------
    % 4) Preallocate component histories and prediction array
    % ------------------------------------------------------------
    tau_hist      = NaN(K,1);   % total trust τ(t)
    tau_lat_hist  = NaN(K,1);   % latent component
    tau_rep_hist  = NaN(K,1);   % reputation component
    tau_sit_hist  = NaN(K,1);   % situational component
    risk_hist     = NaN(K,1);   % risk values
    y_hat         = NaN(nMeas,1); % model predictions at measurement times

    % ------------------------------------------------------------
    % 5) Main simulation loop over time grid
    % ------------------------------------------------------------
    for k = 1:K
        t_k = t_grid(k);

        % Door events scheduled at this grid index (possibly multiple)
        door_indices = doors_at_k{k};

        if ~isempty(door_indices)
            % HRI events: apply each door trial sequentially at time t_k
            for j = 1:numel(door_indices)
                d_idx = door_indices(j);
                D     = doorEvents(d_idx);

                ev = struct();
                ev.t          = t_k;
                ev.type       = "door";
                ev.risk_value = D.risk_value;
                
                switch mode
                    case "simple"
                        ev.outcome    = D.outcome;
                        ev.followed   = D.followed;

                    case "coupled"
                        ev.followed = behavioral_model(state, steepness);
                        if ev.followed ~= D.followed
                            ev.outcome = 1 - D.outcome;
                        else
                            ev.outcome = D.outcome;
                        end
        
                    otherwise
                        error('trust_simulate_or_predict_one_participant: unknown mode "%s". Use "simple" or "coupled".', mode);
                end
                state = trust_step(state, ev, params);
            end
        else
            % No door at this time: advance with pure latent / reputation /
            % situational dynamics.
            ev = struct();
            ev.t    = t_k;
            ev.type = "none";
            state = trust_step(state, ev, params);
        end

        % Store model predictions at any measurement indices assigned to this grid time
        if ~isempty(meas_at_k{k})
            for idx = meas_at_k{k}
                y_hat(idx) = state.tau;
            end
        end

        % Situational trust for logging (using the updated state and last_risk)
        if isnan(state.last_risk)
            tau_sit = 0.0;
        else
            tau_sit = trust_compute_situational( ...
                state.last_risk, state.tau_disp, state.sc, params);
        end

        % Record component histories and risk
        tau_hist(k)     = state.tau;
        tau_lat_hist(k) = state.tau_lat;
        tau_rep_hist(k) = state.tau_rep;
        tau_sit_hist(k) = tau_sit;
        risk_hist(k)    = state.last_risk;
    end

    % ------------------------------------------------------------
    % 6) Package output struct
    % ------------------------------------------------------------
    sim = struct();
    sim.t_grid       = t_grid(:);
    sim.doorEvents   = doorEvents;
    sim.tau_hist     = tau_hist;
    sim.tau_lat_hist = tau_lat_hist;
    sim.tau_rep_hist = tau_rep_hist;
    sim.tau_sit_hist = tau_sit_hist;
    sim.risk_hist    = risk_hist;

    sim.measurements = measurements;
    sim.y_hat        = y_hat;
end
