% By Bo Hu, last edit: 12/28/2024
classdef MultiChoice
    properties
        choice                          % data: choice, should be 1, 2, ..., n_choices
        covariates                      % data: covariates
        covariate_names                 % covariate names
        n_choices                       % number of choice alternatives
        n_covariates                    % number of covariates
        weight                          % sample weights, default is equal weighting
        N                               % sample size
        model_type                      % the model type
        reg_para                        % regularization parameters
        m_gram                          % number of principal components for Gram matrix
        m_hermite                       % degree of Hermite polynomials
        sigma                           % parameter for radial kernel, default is 1/2
        estimated_para                  % the estimated parameters
        pmat                            % fitted choice probabilities
        V                               % pca eigenvectors
        D                               % pca eigenvalues
        covariate_for_pred              % covariate values used for prediction
        pmat_pred                       % predicted probability matrix (from the .predict method)
        hate_mat                        % heterogeneous ATE matrix
    end

    %%%%%%%%%%
    methods
        function obj = MultiChoice(choice_data, covariate_data, ...
                covariate_names, weight)
            % Object constructer 
            %
            % Inputs: 
            % - choice_data: response variable, categories
            % - covariate_data: covariates/covariates/predictors
            % - covariate_names: covariate names
            %   -- covariate names should be in the format of 
            %   -- {'aaa', 'bbb', 'ccc'}. The order of covariate names should
            %   -- be consistent with columns of covariate_data
            % - weight: sample weights. If not provided, use equal weights
            %      as default. 
             
            obj.choice = choice_data;
            obj.covariates = covariate_data;
            obj.covariate_names = covariate_names;
            [obj.N, obj.n_covariates] = size(obj.covariates);
            obj.n_choices = max(max(choice_data));
            % If customized weights are not provided, use equal weights
            if nargin < 4
                weight = ones(obj.N, 1)/obj.N;
            end
            obj.weight = weight;
        end
        
        %%%%%%%%%%
        function obj = estimate(obj, model_type, reg_para, sigma, ...
                ini_guess)
            % Estimation of multi-choice model 
            %
            % Inputs:
            % - model_type: 
            %   -- 'multi_logit': multinomial logit model
            %   -- 'ordered_logit': ordered logit model
            %   -- 'ordered_probit': ordered probit model
            %   -- 'hermite': ordered mulitnomial choice model with random
            %   --      term nonparametrically estimated using Hermite 
            %           polynomials
            % - reg_para: regularization parameter = [m_gram, m_hermite], 
            %           m_hermite is optional, needed only when model_type
            %           is 'hermite'
            %   -- m_gram: number of basis of reproducing kernel Hilbert space
            %   -- m_hermite: degree of Hermite polynomials 
            % - ini_guess: initial guess in optimization algorithm,
            %           optional. 
            %
            % The output trans_data is used in an optional way            
            
            obj.model_type = model_type;
            obj.reg_para = reg_para;

            obj.m_gram = reg_para(1);
            if length(reg_para) == 2
                obj.m_hermite = reg_para(2);
            else
                obj.m_hermite = [];
            end
            
            % kernel trick
            if nargin < 4 || isempty(sigma)
                obj.sigma = 1/2;
            else
                obj.sigma = sigma;
            end
            K = gaussian_kernel(obj.covariates, obj.covariates, obj.sigma);
            [tmp_V, tmp_D] = eigs(K, obj.m_gram);
            % eigenvectors (V) and eigenvalues (D) of Gram matrix
            obj.V = tmp_V;
            obj.D = tmp_D;
            H = obj.V;
            
            % indicator of classes for each individual
            trans_data = zeros(obj.N, obj.n_choices);
            for j = 1:obj.n_choices
                trans_data(:, j) = (obj.choice == j);
            end
            
            % objective function is the -ln L
            objfun = @(x)mlogit_neg_log_likelihood(x, ...
                trans_data, obj.weight, H, obj.model_type, obj.m_hermite);
            
            % initial guess
            % the structure of model parameter is
            % [gram_pc_coeffs; hermite_coeffs; thresholds], in which
            % hermite_coeffs and thresholds are optional, depending on
            % model type. This is a column vector, whose dimension is
            % m_gram + m_hermite + n_choices - 1
            if nargin < 5 
                if strcmp(model_type, "multi_logit")
                    % no hermite coeffs and thresholds
                    ini_guess = zeros(obj.m_gram+1, 2);
                elseif any(strcmp(model_type, {'ordered_probit', ...
                        'ordered_logit'}))
                    % no hermite coeffs 
                    ini_guess = zeros(obj.m_gram+obj.n_choices-1, 1);
                    % the second threshold should be greater than the first one
                    ini_guess(end-obj.n_choices+2:end) = linspace(0, 1, obj.n_choices-1)';
                elseif strcmp(model_type, "hermite")
                    ini_guess = zeros(obj.m_gram+obj.m_hermite+obj.n_choices-1, 1);
                    ini_guess(end-obj.n_choices+2:end) = linspace(0, 1, obj.n_choices-1)';
                else
                    ini_guess = [];
                end
            end
            
            % optimization parameters
            options = optimset('MaxFunEvals',5000000, 'Display', 'off');

            % estimated coeffcients in the kernel expression
            if any(strcmp(model_type, {'multi_logit', 'ordered_probit', ...
                    'ordered_logit'}))
                delta = fminunc(objfun, ini_guess, options);
            elseif strcmp(model_type, 'hermite')
                % constraint function
                confun = @(delta)constraint_function(delta, obj.reg_para);
                delta = fmincon(objfun, ini_guess, [], [], [], [], [], ...
                     [], confun, options);
            end
            obj.estimated_para = delta;

            % fitted transitional probability
            obj.pmat = trans_prob_mat(delta, obj.n_choices, H, ...
                obj.model_type, obj.m_hermite); 
        end

        %%%%%%%%%%
        function obj = predict(obj, covariate_for_pred)
            % Model prediction
            %
            % Inputs: 
            % - covariate_for_pred: values of covariates for prediction

            K_pred = gaussian_kernel(obj.covariates, ...
                covariate_for_pred, obj.sigma);
            H_pred = K_pred * obj.V / obj.D ;
            obj.pmat_pred = trans_prob_mat(obj.estimated_para, obj.n_choices, ...
                H_pred, obj.model_type, obj.m_hermite);
            obj.covariate_for_pred = covariate_for_pred;
             
        end

        %%%%%%%%%% 
        function obj = hate(obj, treatment, contingent, cont_values, ...
                ate_type, contingent_type)
            % Calculate the heterogeneous treatment effect for a single 
            % treatment
            %
            % Inputs:
            % - treatment: treatment variable
            % - contingent: contingent variable
            % - cont_values: contigent variable values at prediction
            %       -- If contingent and cont_values are missing, then 
            %          calculate the average treatment effect
            % - ate_type: 'counterfactual' or 'conditional'
            %       -- If ate_type = 'conditional', need to further 
            %          specify type of contingent variable: 
            %          'continuous' or 'discrete'

            % find location of treament variables
            treat_index = find(ismember(obj.covariate_names, treatment));            
            % create covariate values for treated and untreated group
            treated_group = obj.covariates;
            treated_group(:, treat_index) = ones(obj.N, 1) * ...
                max(obj.covariates(:, treat_index));
            untreated_group = obj.covariates;
            untreated_group(:, treat_index) = ones(obj.N, 1) * ...
                min(obj.covariates(:, treat_index));
            
            % if contingent information is not given, calculate ATE.
            % Otherwise, calculate CATE
            if nargin > 2
                cont_index = find(ismember(obj.covariate_names, ...
                    contingent));
                n_eval = length(cont_values);
            else
                n_eval = 1;
            end
            
            % calculate CATE
            obj.hate_mat = zeros(n_eval, obj.n_choices);
            for i = 1:n_eval
                if nargin > 2
                    treated_group(:, cont_index) = ones(obj.N, 1) * ...
                        cont_values(i);
                    untreated_group(:, cont_index) = ones(obj.N, 1) * ...
                        cont_values(i);
                end
                pmat_1 = obj.predict(treated_group).pmat_pred;
                pmat_0 = obj.predict(untreated_group).pmat_pred;
                if strcmp(ate_type, 'counterfactual')
                    obj.hate_mat(i, :) = obj.weight' / sum(obj.weight) ...
                        * (pmat_1 - pmat_0);
                elseif strcmp(ate_type, 'conditional')
                    if strcmp(contingent_type, 'continuous')
                        % If contingent variable is continuous, use kernel
                        % method to estimate conditional expectation with 
                        % Gaussian kernel
                        kernel_weights = exp(- (obj.covariates(:, ...
                            cont_index) - cont_values(i)).^2 / 2);
                        obj.hate_mat(i, :) = (obj.weight.*kernel_weights)' ...
                            * (pmat_1 - pmat_0) / (sum(obj.weight.*kernel_weights));
                    end
                end
            end
        end

        %%% Treatment and Heterogeneous Treatment Effect (Multiple treatment)
        function obj = hate_multi_treat(obj, treatment, treat_values, ...
                control_values, contingent, cont_values, ...
                ate_type, contingent_type, class_thresholds)
            % Calculate the heterogeneous treatment effect for multiple 
            % treatments. Need to specify treatment values
            % 
            % Inputs:
            % - treatment: treatment variable
            % - treat_values: treatment variable values for treated group,
            %       should be row vector
            % - control_values: treatment variable values for control group,
            %       should be row vector
            % - contingent: contingent variable
            % - cont_values: contigent variable values at prediction
            %       -- If contingent and cont_values are missing, then 
            %          calculate the average treatment effect            
            % - ate_type: 'counterfactual' or 'conditional'
            %       -- If ate_type = 'conditional', need to further 
            %          specify type of contingent variable: 'continuous' 
            %          or 'continuous-class' or 'discrete'. If 
            %          'continuous-class', have to provide class
            %          thresholds as a m-dimensional vector, then there 
            %          are m+1 class

            % find location of treament and contingent variables
            treat_index = zeros(length(treatment), 1);
            for i = 1:length(treatment)
                treat_index(i) = find(ismember(obj.covariate_names, treatment(i)));
            end      
            treated_group = obj.covariates;
            treated_group(:, treat_index) = ones(obj.N, 1) * treat_values;
            untreated_group = obj.covariates;
            untreated_group(:, treat_index) = ones(obj.N, 1) * ...
                control_values;
               
            if nargin > 4
                cont_index = find(ismember(obj.covariate_names, contingent));
                n_eval = length(cont_values);
            else
                n_eval = 1;
            end
            
            obj.hate_mat = zeros(n_eval, obj.n_choices);
            for i = 1:n_eval
                if nargin > 4 && (strcmp(ate_type, 'counterfactual') ...
                        || strcmp(contingent_type, 'continuous'))
                    % use given contingent values only when contingent is
                    % continuous. Otherwise, use data contingent values
                    % since 
                    treated_group(:, cont_index) = ones(obj.N, 1) * ...
                        cont_values(i);
                    untreated_group(:, cont_index) = ones(obj.N, 1) * ...
                        cont_values(i);
                end
                pmat_1 = obj.predict(treated_group).pmat_pred;
                pmat_0 = obj.predict(untreated_group).pmat_pred;
                if strcmp(ate_type, 'counterfactual')
                    obj.hate_mat(i, :) = obj.weight' / sum(obj.weight) ...
                        * (pmat_1 - pmat_0);
                elseif strcmp(ate_type, 'conditional')
                    if strcmp(contingent_type, 'continuous')
                        kernel_weights = exp(- (obj.covariates(:, cont_index) ...
                            - cont_values(i)).^2 / 2);
                        obj.hate_mat(i, :) = (obj.weight.*kernel_weights)' ...
                            * (pmat_1 - pmat_0) / (sum(obj.weight.*kernel_weights));
                    elseif strcmp(contingent_type, 'continuous-class')
                        num_classes = length(class_thresholds) + 1;
                        obj.hate_mat = zeros(num_classes, self.n_choices);
                        pmat_diff = pmat_1 - pmat_0;
                        % the first class
                        idx_left = (obj.covariates(:, cont_index) <= ...
                            class_thresholds(1));
                        pmat_diff_left = pmat_diff(idx_left, :);
                        weight_left = obj.weight(idx_left);
                        obj.hate_mat(1, :) = weight_left' / ...
                            sum(weight_left) * pmat_diff_left;
                        % the middle classes
                        if num_classes > 2
                            for k = 2:(num_classes-1)
                                idx_middle = ((obj.covariates(:, cont_index) > ...
                                    class_thresholds(k-1)) ...
                                    & (obj.covariates(:, cont_index) <= ...
                                    class_thresholds(k)));
                                pmat_diff_middle = pmat_diff(idx_middle, :);
                                weight_middle = obj.weight(idx_middle);
                                obj.hate_mat(k, :) = weight_middle' / ...
                                    sum(weight_middle) * pmat_diff_middle;
                            end
                        end
                        % the last class
                        idx_right = (obj.covariates(:, cont_index) > ...
                            class_thresholds(end));
                        pmat_diff_right = pmat_diff(idx_right, :);
                        weight_right = obj.weight(idx_right);
                        obj.hate_mat(end, :) = weight_right' / ...
                            sum(weight_right) * pmat_diff_right;
                    end
                end
            end
        end

        %%%%%%%%%%
        function obj = loo_cross_validation(obj, model_type_cv)
            % Conduct leave-one-out cross validation to determine tuning
            % parameters
            %
            % Inputs:
            % - model_type_cv: the model type
            %

            % maximum m_gram and m_hermite
            m_gram_try = 8;
            m_hermite_try = 3;
            loss_mse_weight = zeros(m_gram_try, m_hermite_try);

            for m_cv = 1:m_gram_try
                for J_cv = 1:m_hermite_try
                    cum_loss_mse_weight = 0;
                    for i = 1:obj.N            
                        fprintf('m_cv = %d, J_cv = %d, i = %d \n', ...
                            m_cv, J_cv, i);

                        choice_cv = obj.choice;
                        covariates_cv = obj.covariates;
                        weight_cv = obj.weight;
                        % leave one out
                        choice_cv(i, :) = [];
                        covariates_cv(i, :) = [];
                        weight_cv(i, :) = [];
            
                        % estimation and prediction
                        reg_para_cv = [m_cv; J_cv];
                        pred_dist = MultiChoice(choice_cv, covariates_cv, ...
                            obj.covariate_names, weight_cv).estimate(...
                            model_type_cv, reg_para_cv).predict( ...
                            obj.covariates(i, :)).pmat_pred;

                        data_dist = zeros(1, self.n_choices);
                        data_dist(obj.choice(i)) = 1;
        
                        % with weights
                        cum_loss_mse_weight = cum_loss_mse_weight + ...
                            obj.weight(i) * sum((data_dist - pred_dist).^2);
                    end
        
                    loss_mse_weight(m_cv, J_cv) = cum_loss_mse_weight/...
                        sum(obj.weight);
                end
            end
            
            % find the positions of the 5 smallest numbers
            k = 5;
            % Flatten the matrix and sort it to find the positions of k 
            % smallest elements
            [~, sortedIndices] = sort(loss_mse_weight(:));
            smallestKIndices = sortedIndices(1:k);
            
            % Convert linear indices to row and column indices
            [rowIndices, colIndices] = ind2sub(size(loss_mse_weight), ...
                smallestKIndices);
            
            % Display the positions of the k smallest numbers
            fprintf('Positions of the %d smallest numbers:\n', k);
            for i = 1:k
                fprintf('LOOCV selects (m, r) = (%d, %d)\n', ...
                    rowIndices(i), colIndices(i));
            end
        end
    end
end
