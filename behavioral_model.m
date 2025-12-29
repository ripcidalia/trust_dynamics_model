function behavior = behavioral_model(state, k)
% behavioral_model  Predict binary human action given current trust state.
%
%   behavior = behavioral_model(state)
%
% This function runs an instance of the probabilistic behavioral model
% for a given trust state and steepness k of the probabilistic transition.
% It models the human action (follow/not follow) as a Bernoulli random
% variable with probability p. The decision boundary is modeled as a smooth
% continuous, sigmoid function of the difference between the trust state
% and the decision threshold:
%
%       p = 1/(1 + exp(-k*(tau-self_confidence))).
%
% Inputs:
%   state     Current trust state with (relevant) fields:
%               .tau         Current total trust in [0,1]
%               .sc          Self-confidence (scalar in [0,1])
%
%   k         Controls the steepness of the probabilistic transition between the two
%             behavioral actions (follow/ not follow). If omitted or empty, steepness
%             defaults to 10.
%
% Output:
%   behavior  Predicted binary action: 1 = followed, 0 = not followed

    if nargin < 2 || isempty(k)
        k = 10;
    end

    if ~isfield(state, "tau") || isempty(state.tau)
        error("behavioral_model: no trust value specified.");
    else
        tau = state.tau;
    end

    if ~isfield(state, "sc") || isempty(state.sc)
        warning("behavioral_model: no self-confidence value specified. Proceeding with default 0.5");
        self_confidence = 0.5;
    else
        self_confidence = state.sc;
    end

    exp_term = -k*(tau - self_confidence);

    p = 1/(1 + exp(exp_term));

    if rand() < p
        behavior = 1;
    else
        behavior = 0;
    end

