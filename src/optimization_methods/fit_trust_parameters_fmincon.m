function [theta_hat, fval, exitflag, output] = fit_trust_parameters_fmincon(dt, theta0)
% fit_trust_parameters_fmincon  Fit global trust-dynamics parameters via WLS + fmincon.
%
%   [theta_hat, fval, exitflag, output] = ...
%       fit_trust_parameters_fmincon(dt)
%
%   [theta_hat, fval, exitflag, output] = ...
%       fit_trust_parameters_fmincon(dt, theta0)
%
% Inputs:
%   dt     - time step for simulation (seconds)
%   theta0 - optional initial parameter guess (7x1 vector). If not provided,
%            a default heuristic initialisation is used.
%
% Outputs:
%   theta_hat - estimated parameter vector:
%       [ lambda_rep
%         phi_fail
%         psi_succ
%         a_succ
%         lambda_sit
%         lambda10
%         kappa01
%         theta_sit ]
%
%   fval      - final WLS cost
%   exitflag  - fmincon exit flag
%   output    - fmincon output struct

    if nargin < 1 || isempty(dt)
        error('fit_trust_parameters_fmincon requires dt.');
    end

    %% ------------------------------------------------------------
    % 1) Load participants and measurement weights
    % -------------------------------------------------------------
    fprintf('[fit_trust_parameters_fmincon] Loading data...\n');

    matP = "derived/participants_probes_mapped_stepM4.mat";
    if ~isfile(matP)
        error("File not found: %s. Run preprocessing (M1–M4) first.", matP);
    end
    S = load(matP, "participants_probes_mapped");
    if ~isfield(S, "participants_probes_mapped")
        error("Variable participants_probes_mapped not found in %s.", matP);
    end
    participants = S.participants_probes_mapped;
    N = numel(participants);
    fprintf('[fit_trust_parameters_fmincon] Loaded %d participants.\n', N);

    % Measurement weights (already computed via LOPO CV)
    matW = "derived/measurement_weights.mat";
    if ~isfile(matW)
        error("File not found: %s. You said this should contain 'weights'.", matW);
    end
    W = load(matW, "weights");
    if ~isfield(W, "weights")
        error("measurement_weights.mat does not contain 'weights'.");
    end
    meas_weights = W.weights;

    fprintf('[fit_trust_parameters_fmincon] Loaded measurement weights:\n');
    disp(meas_weights);

    %% ------------------------------------------------------------
    % 2) Initial parameter guess theta0
    % -------------------------------------------------------------
    % Order:
    %   1: lambda_rep
    %   2: phi_fail
    %   3: psi_succ
    %   4: a_succ
    %   5: lambda_sit
    %   6: lambda10
    %   7: kappa01
    %   8: theta_sit

    if nargin < 2 || isempty(theta0)
        theta0 = [
            1e-3;   % lambda_rep
            0.15;   % phi_fail
            0.10;   % psi_succ
            0.20;   % a_succ
            1.0;    % lambda_sit
            1e-4;   % lambda10
            1e-4;   % kappa01
            0.5;    % theta_sit
        ];
        fprintf('[fit_trust_parameters_fmincon] Using default initial guess.\n');
    else
        fprintf('[fit_trust_parameters_fmincon] Using user-provided initial guess.\n');
        theta0 = theta0(:);
    end

    %% ------------------------------------------------------------
    % 3) Bounds on parameters 
    % -------------------------------------------------------------
    %   lambda_rep ∈ [1e-4, 0.1]
    %   phi_fail   ∈ [0.01, 0.99]
    %   psi_succ   ∈ [0.01, 0.99], with constraint psi_succ <= phi_fail
    %   a_succ     ∈ [0.01, 0.65]
    %   lambda_sit ∈ [0.1, 2]
    %   lambda10   ∈ [1e-6, 1e-2]
    %   kappa01    ∈ [1e-6, 1e-2]
    %   theta_sit  ∈ [0, 1]

    lb = [
        1e-4;    % lambda_rep
        0.01;    % phi_fail
        0.01;    % psi_succ
        0.01;    % a_succ
        0.10;    % lambda_sit
        1e-6;    % lambda10
        1e-6;    % kappa01
        0;       % theta_sit
    ];

    ub = [
        0.10;    % lambda_rep
        0.99;    % phi_fail
        0.99;    % psi_succ
        0.65;    % a_succ
        2.0;     % lambda_sit
        1e-2;    % lambda10
        1e-2;    % kappa01
        1;       % theta_sit
    ];

    %% ------------------------------------------------------------
    % 4) Define SAFE objective wrapper for fmincon
    % -------------------------------------------------------------
    % Raw cost:
    raw_fun = @(theta) trust_cost_all(theta, participants, dt, meas_weights);

    % Robust wrapper: catches errors AND non-real/non-scalar/NaN/Inf and
    % returns a large penalty so fmincon never sees an undefined objective.
    function J = safe_obj(theta_vec)
        persistent n_error_calls
        if isempty(n_error_calls)
            n_error_calls = 0;
        end

        try
            J_raw = raw_fun(theta_vec);

            % ---- Hard checks: J must be a real, finite scalar ----
            if ~isscalar(J_raw) || ~isreal(J_raw) || ~isfinite(J_raw)
                n_error_calls = n_error_calls + 1;
                if n_error_calls <= 10
                    warning(['safe_obj: invalid cost value (non-scalar, ' ...
                             'non-real, or non-finite). Using large penalty.']);
                end
                J = 1e6;
                return;
            end

            % All good
            J = J_raw;

        catch ME
            n_error_calls = n_error_calls + 1;
            if n_error_calls <= 10
                warning('safe_obj: error in trust_cost_all at some theta: %s', ME.message);
            end
            J = 1e6;
        end
    end

    obj = @(theta) safe_obj(theta);

    %% ------------------------------------------------------------
    % 5) Nonlinear constraints (phi_fail >= psi_succ)
    % -------------------------------------------------------------
    % c(theta) <= 0
    %   -> c1 = psi_succ - phi_fail <= 0
    %
    % No equality constraints.
    %
    function [c, ceq] = trust_nonlcon(theta_vec)
        phi_fail = theta_vec(2);
        psi_succ = theta_vec(3);

        % Enforce psi_succ <= phi_fail
        c1 = psi_succ - phi_fail;

        c   = c1;     % inequality vector
        ceq = [];     % no equalities
    end

    nonlcon = @trust_nonlcon;

    %% ------------------------------------------------------------
    % 6) Check initial cost at theta0
    % -------------------------------------------------------------
    J0_raw = raw_fun(theta0);
    fprintf('[fit_trust_parameters_fmincon] Initial raw cost at theta0: %g\n', J0_raw);

    if ~isfinite(J0_raw)
        warning(['Initial theta0 gives non-finite cost (%g). ' ...
                 'safe_obj will use a large penalty instead.'], J0_raw);
    end

    %% ------------------------------------------------------------
    % 7) Call fmincon
    % -------------------------------------------------------------
    options = optimoptions('fmincon', ...
        'Algorithm', 'interior-point', ...
        'Display', 'iter', ...
        'MaxIterations', 200, ...
        'MaxFunctionEvaluations', 2000, ...
        'StepTolerance', 1e-8, ...
        'OptimalityTolerance', 1e-6);

    fprintf('[fit_trust_parameters_fmincon] Starting fmincon...\n');

    [theta_hat, fval, exitflag, output] = fmincon( ...
        obj, theta0, [], [], [], [], lb, ub, nonlcon, options);

    fprintf('[fit_trust_parameters_fmincon] Done.\n');
    fprintf('  exitflag = %d\n', exitflag);
    fprintf('  final cost fval = %.6f\n', fval);
    fprintf('  theta_hat =\n');
    disp(theta_hat);

end
