function log_like = mlogit_neg_log_likelihood(delta, data, weight, H, model_type, m_hermite)
% Inputs:
%  - delta: model parameters
%  - data: data (transition indicators)
%  - weight: sample weights
%  - H: the H matrix from Gram matrix
%  - model_type: type of model, see MultiChoice()
%  - m_hermite: degree of Hermite polynomial for the random term
%       If model_type = 'hermite', m_hermite is necessary, other wise just supply []

% obtain the transition prbability matrix
p_mat = trans_prob_mat(delta, size(data, 2), H, model_type, m_hermite);

% negative log likelihood, to accommodate with minimization routine
log_like = -sum(sum(weight .* data .* log(p_mat)));

end