function [c, ceq] = constraint_function(delta, reg_para)
% For hermite polynomial approach, should supply reg_para, whose second 
% component is the order of Hermite polynomials used

m_gram = reg_para(1);
m_hermite = reg_para(2);

tau = delta(m_gram+1:m_gram+m_hermite);

gam = GNaux_gamma(tau);
const_a = GNaux_a(tau, 1);

c = [];
ceq = gam' * const_a(2:end);

