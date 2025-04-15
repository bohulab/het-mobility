function const_a = GNaux_a(tau, flag_extra_term)
% auxilary functions: genenrate the a constants
% input:
%   tau: parameter for Gallant & Nychka form. A column vector

R = length(tau);

if flag_extra_term == 1
    length_a = 2 * R + 2;
else
    length_a = 2 * R + 1;
end

const_a = zeros(length_a,1);

const_a(1) = 1;
const_a(2) = 0;
const_a(3) = 1;

for h=3:(length_a-1)
    const_a(h+1) = (h-1)*const_a(h-1);
end

end