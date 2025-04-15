function F = GNdist(x, tau)
% distribution function of Gallant & Nychka (1987) form
% input:
%   x: evaluation point
%   tau: parameter for Gallant & Nychka form. A column vector
%   

% augment tau so that the first element of tau is 1 (normalization)
n = length(x);
R = length(tau);

gam = GNaux_gamma(tau);

% calculate the a constants
const_a = GNaux_a(tau, 0);

% calculate the A functions evaluated at x
func_A = zeros(n, 2*R+1);
func_A(:, 1) = normcdf(x);
func_A(:, 2) = -normpdf(x);
func_A(:, 3) = x .* func_A(:, 2) + func_A(:, 1);

for h=3:2*R
    func_A(:, h+1) = x .* (func_A(:, h) -(h-2)*func_A(:, h-2) )  + (h-1)*func_A(:, h-1);
end

% calculate psi
psi = gam' * const_a;

% calculate distribution function F
F = func_A * gam / psi;

end