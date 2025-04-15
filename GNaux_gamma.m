function gam = GNaux_gamma(tau)
% auxilary functions: genenrate the gamma's 
% input:
%   tau: parameter for Gallant & Nychka form. A column vector

R = length(tau);
ctau = [1; tau];

% calculate gamma's (borrowed from Guo. need to confirm)
T = ( ctau(end:-1:1) )*ctau';
gam = zeros(2*R+1,1);

for h=0:2*R
    gam(h+1) = sum( diag(T,h-R) );
end

end