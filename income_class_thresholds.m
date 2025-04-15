function [upper, lower] = income_class_thresholds(income, weight)
% calculate the threshold for income classes according to Pew Research 
% Center's criteria

median_val = weightedMedian(income, weight);
upper = median_val * 2;
lower = median_val * 2/3;

end