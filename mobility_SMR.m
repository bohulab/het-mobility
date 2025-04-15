% mobility 
% By Bo Hu, last edit: 12/28/2024

%% load data
clear
close all
rng(2024)       % Fix random generator seed for replicability
conduct_bootstrap = 0;      % Bootstrap can take time
plot_results = 0;

%% read data
data = readmatrix("child_parent_data.xls");
% For robustness check
% data = readmatrix("child_parent_data_10_15.xls");
[N, ~] = size(data);
est_weight = data(:, end);

% Child income class
% Use Pew Research Center's categorization
chdinc_classes_obj = Income(data(:, 1), est_weight).classify();
% Use quantiles for categorization
% chdinc_classes_obj = Income(data(:, 1), est_weight).classify_quantiles(5);
y_data = chdinc_classes_obj.income_class;
n_choices = max(y_data);

% Covariates
covs = data(:, 2:end-1);
covs(:, 1) = log(covs(:, 1));
% Normalization in order to construct Gram matrix
[x_data, cov_mean, cov_std] = data_normalize(covs);
factor_names = {'prt_inc', 'sex', 'college', 'race', 'prt_age', ...
    'birthweight', 'prt_college'};

%% model estimation
% Tuning parameters
m_gram = 7;
m_hermite = 2;
reg_para = [m_gram; m_hermite];
model_type = 'hermite';
% Other possible model types for robustness
% model_type = 'ordered_probit';
% model_type = 'ordered_logit';

% Estimate the model
mb_model = MultiChoice(y_data, x_data, factor_names, est_weight)...
    .estimate(model_type, reg_para);

%% Treatment effects
% Get the contingent variable evaluation points
grid_income_pred = (6:0.1:12)';         % parental income
grid_prtage_pred = (18:0.2:45)';        % parental age at child's birth
grid_income_pred_norm = (grid_income_pred - cov_mean(1))/cov_std(1);
grid_prtage_pred_norm = (grid_prtage_pred - cov_mean(5))/cov_std(5);

% Single treatment effects
hate_sex = mb_model.hate('sex', 'prt_inc', ...
    grid_income_pred_norm, 'conditional', 'continuous').hate_mat;
hate_prtedu = mb_model.hate('prt_college', 'prt_inc', ...
    grid_income_pred_norm, 'conditional', 'continuous').hate_mat;
hate_race = mb_model.hate('race', 'prt_inc', ...
    grid_income_pred_norm, 'conditional', 'continuous').hate_mat;
hate_prtedu_age = mb_model.hate('prt_college', 'prt_age', ...
    grid_prtage_pred_norm, 'conditional', 'continuous').hate_mat;

% Multiple treatment effect
treat_vars = {'race', 'prt_college', 'prt_age'};
treat_position = zeros(length(treat_vars), 1);
for i = 1:length(treat_vars)
    treat_position(i) = find(ismember(factor_names, treat_vars(i)));
end

treat_vals = zeros(1, 3);
control_vals = zeros(1, 3);

% Treat and control values for race 
% Advantaged group: black = 0
% Disadvantaged group: black = 1
treat_vals(1) = min(x_data(:, treat_position(1)));
control_vals(1) = max(x_data(:, treat_position(1)));

% Treat and control values for prt_college
% Advantaged group: prt_college = 1
% Disadvantaged group: prt_college = 0
treat_vals(2) = max(x_data(:, treat_position(2)));
control_vals(2) = min(x_data(:, treat_position(2)));

% Treat and control values for prt_age
% Advantaged group: prt_age = 30
% Disadvantaged group: prt_age = 18
treat_vals(3) = (30-cov_mean(treat_position(3)))/cov_std(treat_position(3));
control_vals(3) = (18-cov_mean(treat_position(3)))/cov_std(treat_position(3));

hate_multiple = mb_model.hate_multi_treat(treat_vars, ...
    treat_vals, control_vals, 'prt_inc', grid_income_pred_norm, ...
    'conditional', 'continuous').hate_mat;

%% bootstrap for confidence bands
B = 1000;
n_pred_income = length(grid_income_pred);
n_pred_prtage = length(grid_prtage_pred);

hate_sex_boot = nan(n_pred_income, n_choices, B);
hate_race_boot = nan(n_pred_income, n_choices, B);
hate_prtedu_boot = nan(n_pred_income, n_choices, B);
hate_prtedu_age_boot = nan(n_pred_prtage, n_choices, B);
hate_multiple_boot = nan(n_pred_income, n_choices, B);

if conduct_bootstrap == 1
    % Set initial guess
    delta_0 = mb_model.estimated_para;
    
    for b = 1:B
        fprintf('b = %d \n', b)
        % Resample
        boot_ind = randsample(N, N, true, ...
            est_weight/sum(est_weight));
        model_boot = MultiChoice(y_data(boot_ind), x_data(boot_ind, :), ...
            factor_names).estimate(model_type, reg_para, [], delta_0);
        
        % Calculate treatment effect
        hate_sex_boot(:, :, b) = model_boot.hate('sex', 'prt_inc', ...
            grid_income_pred_norm, 'conditional', 'continuous').hate_mat;
        hate_race_boot(:, :, b) = model_boot.hate('race', 'prt_inc', ...
            grid_income_pred_norm, 'conditional', 'continuous').hate_mat;
        hate_prtedu_boot(:, :, b) = model_boot.hate('prt_college', 'prt_inc', ...
            grid_income_pred_norm, 'conditional', 'continuous').hate_mat;
        hate_prtedu_age_boot(:, :, b) = model_boot.hate('prt_college', ...
            'prt_age', grid_prtage_pred_norm, 'conditional', ...
            'continuous').hate_mat;
        hate_multiple_boot(:, :, b) = model_boot.hate_multi_treat(treat_vars, ...
            treat_vals, control_vals, 'prt_inc', grid_income_pred_norm, ...
            'conditional', 'continuous').hate_mat;
    end

    save('./results/hate_boot_results.mat', ...
        'hate_sex_boot', 'hate_race_boot', 'hate_prtedu_boot', ...
        'hate_prtedu_age_boot', 'hate_multiple_boot');
end

%% plot the results
if plot_results == 1

    load('./results/hate_boot_results.mat')

    hate.hate_sex = hate_sex;
    hate.hate_race = hate_race;
    hate.hate_prtedu = hate_prtedu;
    hate.hate_prtedu_age = hate_prtedu_age;
    hate.hate_multiple = hate_multiple;
    
    hate_boot.hate_sex_boot = hate_sex_boot;
    hate_boot.hate_race_boot = hate_race_boot;
    hate_boot.hate_prtedu_boot = hate_prtedu_boot;
    hate_boot.hate_prtedu_age_boot = hate_prtedu_age_boot;
    hate_boot.hate_multiple_boot = hate_multiple_boot;
    
    report_vars = {'sex', 'race', 'prtedu', 'prtedu_age', 'multiple'};
    x_labels = {'parental income', 'parental income', ...
        'parental income', 'parental age at birth', 'parental income'};
    
    for i = 1:length(report_vars)
        plot_name = ['prob_diff_' report_vars{i}];
        est_y_var = ['hate_' report_vars{i}];
        boot_y_var = ['hate_' report_vars{i} '_boot'];
    
        if strcmp(x_labels{i}, 'parental income')
            x_vals = grid_income_pred;
        elseif strcmp(x_labels{i}, 'parental age at birth')
            x_vals = grid_prtage_pred;
        end
     
        ate_plot(plot_name, x_labels{i}, x_vals, ...
            hate.(est_y_var), hate_boot.(boot_y_var), 'q5'); 
    end
end