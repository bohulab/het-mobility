% By Bo Hu, last edit: 12/28/2024
classdef Income
    % the Income class
    properties
        N                       % data sample size
        income_data             % income data
        weight                  % data weight
        income_class            % income class (1=low, 2=middle, 3=high) 
        income_class_dist       % income class distribution
        lower                   % upper bound for low income class
        upper                   % lower bound for high income class
        quantiles               % quantiles of income distribution (evenly 
                                %   spaced quantiles, starts after the 0 
                                %   quantile and ends before the 1 quantile)
    end

    methods
        %%%
        function obj = Income(income_data, weight)
            % Class constructor
            % The input income_data should be a column vector
            [obj.N, m] = size(income_data);
            
            % Conduct some consistency check
            
            if m ~= 1
                error('Income have multiple columns.')
            elseif nargin == 2
                [w_m, w_n] = size(weight);
                if w_n ~= 1
                    error('Weight have multiple columns.')
                elseif obj.N ~=w_m
                    error('Dimensions of income and weight do not match.')
                end
            end

            % Weight is optional. If missing, assign equal weights
            if nargin == 1
                weight = ones(obj.N, 1)/obj.N;
            end

            obj.income_data = income_data;
            obj.weight = weight;
        end

        %%%
        function obj = classify(obj, lower, upper)
            % classify incomes into low-, middle-, and high-classes
            if nargin < 2
                [upper, lower] = income_class_thresholds(obj.income_data, obj.weight);
            end
            income_rank = 2 * ones(obj.N, 1);
            income_rank(obj.income_data < lower) = 1;
            income_rank(obj.income_data > upper) = 3;
            obj.income_class = income_rank;
            obj.lower = lower;
            obj.upper = upper;
        end

        %%%
        function obj = classify_quantiles(obj, n_groups)
            % classify incomes into groups according to quantiles
            % e.g, n_groups = 4 means we classify income into 4 groups
            %       based on quartiles               
            
            % get quantiles
            qts = quantile(obj.income_data, n_groups-1);
            % set highest rank
            income_rank = n_groups * ones(obj.N, 1);
            for i = n_groups-1:-1:1
                income_rank(obj.income_data < qts(i)) = i;
            end
            obj.income_class = income_rank;
            obj.quantiles = qts;

        end

        %%%
        function obj = class_distribution(obj)
            % calculate the distribution of income classes
            % obtain the number of income groups
            n_groups = max(obj.income_class);
            inc_dist = zeros(n_groups, 1);
            for i = 1:n_groups
                inc_dist(i) = sum(obj.weight(obj.income_class == i))...
                    / sum(obj.weight);             
            end
            obj.income_class_dist = inc_dist; 
        end

    end
end