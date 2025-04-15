function [normed_data, data_mean, data_std] = data_normalize(data)
% function used to normalize data    
data_mean = mean(data);
data_std = std(data);
normed_data = (data - data_mean) ./ data_std;
end