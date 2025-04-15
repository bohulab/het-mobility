function p_mat = trans_prob_mat(delta, n_choices, H, model_type, m_hermite)
% Obtain the transition probability matrix
% Input
%  - delta: model parameters
%  - data: data (transition indicators)
%  - H: the H matrix from Gram matrix
%  - model_type: type of model, see MultiChoice()
%  - m_hermite: degree of Hermite polynomial for the random term
%       If model_type = 'hermite', m_hermite is necessary, other wise just supply []
% output
%   p_mat: 1*n_choices vector, n_choices is the number of columns of data
%
% By Bo Hu, last edit: 12/28/2024

% sample size
[N, ~] = size(H);

% The case of multivariate logit model
if strcmp(model_type, "multi_logit")
    % The multivariate logit model sets the probability of falling into 
    % category j to be e^{g_j(x)}/ (\sum_k e^{g_k(x)}), j= 1, ..., J
    %
    % In applications, if we have J categories, we should have J
    % g_j functions corresponding to the J categories. But for
    % identification in the multivariate logit model, the first function is
    % set to be the constant function zero. Therefore we only need to 
    % estimate the remaining m-1 g_j's. That is, delta has J-1 columns
    %
    % We include a constant term in delta
    exp_g = exp(H * delta(1:end-1, :) + delta(end, :));

    p_mat_num = [ones(N, 1) exp_g]; % probability numerator
    p_mat_den = sum(p_mat_num, 2) * ones(1, n_choices);
    p_mat = p_mat_num./p_mat_den;

elseif strcmp(model_type, "ordered_probit")
    % The ordered logit model sets the probability of falling into category
    % j to be F(\gamma_j - g(x))-F(\gamma_{j-1} - g(x))), j = 1, ..., J-1
    % with \gamma_0 = -infty, \gamma_J = infty, and F is the distribution
    % function of the standard normal distribution
    %
    % For identification, we should not include constant term in g, or we
    % include a constant term in g, but also fix one of the \gamma_j for
    % j=1, 2, ..., J-1.
    % 
    % delta therefore contains m_gram kernel coefficients and n_choices-1 
    % thresholds, which is set as a column vector
    p_mat = zeros(N, n_choices);
    p_cum = zeros(N, 1);
    for i = 1:n_choices-1
        p_mat(:, i) = normcdf(delta(end+1-n_choices+i, :) ...
            - H*delta(1:end-n_choices+1, :)) - p_cum;
        p_cum = p_cum + p_mat(:, i);
    end
    p_mat(:, n_choices) = 1 - p_cum;

elseif strcmp(model_type, "ordered_logit")
    % exactly the same as the ordered_probit, but distribution function
    % replaced by the logit distribution function
    p_mat = zeros(N, n_choices);
    p_cum = zeros(N, 1);
    for i = 1:n_choices-1
        p_mat(:, i) = logitcdf(delta(end+1-n_choices+i, :) ...
            - H*delta(1:end-n_choices+1, :)) - p_cum;
        p_cum = p_cum + p_mat(:, i);
    end
    p_mat(:, n_choices) = 1 - p_cum;

elseif strcmp(model_type, "hermite")
    % Hermite polynomial approach as in Gallant and Nychka (1987)
    p_mat = zeros(N, n_choices);
    p_cum = zeros(N, 1);
    for i = 1:n_choices-1
        p_mat(:, i) = GNdist(delta(end+1-n_choices+i, :) ...
            - H*delta(1:end-m_hermite-n_choices+1, :), ...
            delta(end-m_hermite-n_choices+2:end-n_choices+1, :)) - p_cum;
        p_cum = p_cum + p_mat(:, i);
    end
    p_mat(:, n_choices) = 1 - p_cum;
end

end