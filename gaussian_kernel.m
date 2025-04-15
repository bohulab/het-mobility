function K = gaussian_kernel(x, y, kappa)
% Output the Gram matrix
%
% Expression of Gaussian kernel: k(t, s) = exp(-kappa*|t-s|^2)
% Matrix expression with data K = [k_ij] = [k(y_i, x_j)]
% 
% Input:
% - x: data point where kernel basis is formed
% - y: points at which kernel basis are evaluated
%       Both x and y: each row vector corresponds to one observation
%       row vectors of x and y are x_j' and y_i'

K = exp(-kappa*(diag(y*y') + diag(x*x')' - 2*y*x'));

end