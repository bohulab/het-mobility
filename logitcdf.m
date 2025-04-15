function y = logitcdf(x)
% CDF function for logit 
% x could be scalar or array
    y = exp(x)./(1+exp(x));
end