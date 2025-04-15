function ate_plot(plot_name, x_lab, x, est_y, boot_y, extra_label)
% Input:
% - plot_name: plot name (string)
% - x_lab = x axis label (string)
% - x: column vector
% - est_y: estimated y vector, each line in a column
% - y_boot: bootstrapped y vector
% - para_labels: extra lables for naming plots,
%   - use 'q3', 'q4', 'q5' to when threshold is set as quantiles

boot_y_mean = mean(boot_y, 3);
boot_y_025 = quantile(boot_y, 0.025, 3) - boot_y_mean;
boot_y_050 = quantile(boot_y, 0.05, 3) - boot_y_mean;
boot_y_950 = quantile(boot_y, 0.95, 3) - boot_y_mean;
boot_y_975 = quantile(boot_y, 0.975, 3) - boot_y_mean;
% plots parameters
% xlim parameter
xlim_para = [min(x) max(x)];
y_max = max(max(est_y + boot_y_975));
y_min = min(min(est_y + boot_y_025));
ylim_para = [y_min*1.1 y_max*1.3];

% get the number of income classes
n_class = size(boot_y_mean, 2);
plot_rows = ceil((n_class+1)/2);

% subplot titles
if strcmp(extra_label, 'q3')
    title_str = {'1st Tertile', '2nd Tertile', '3rd Tertile'};
elseif strcmp(extra_label, 'q4')
    title_str = {'1st Quartile', '2nd Quartile', '3rd Quartile', '4th Quartile'};
elseif strcmp(extra_label, 'q5')
    title_str = {'1st Quintile', '2nd Quintile', '3rd Quintile', '4th Quintile', '5th Quintile'};
else
    title_str = {'Low', 'Middle', 'High'};
end

line_color = [[224 143 51]/255;
    [246 192 66]/255;
    [170 179 58]/255;
    [145, 172, 223]/255;
    [159, 124, 196]/255];

% plots
figure('Name', plot_name)
% plot effects to each class on one plot
subplot(plot_rows, 2, 1)
plot(x, est_y(:, 1), 'LineWidth', 2, 'Color', line_color(1, :));
hold on;
plot(x, est_y(:, 2), 'LineWidth', 2, 'Color', line_color(2, :), 'LineStyle', '--');
plot(x, est_y(:, 3), 'LineWidth', 2, 'Color', line_color(3, :), 'LineStyle', '-.');
if n_class >= 4
    plot(x, est_y(:, 4), 'LineWidth', 2, 'Color',line_color(4, :), 'LineStyle',':')
end
if n_class >=5
    plot(x, est_y(:, 5), 'LineWidth', 2, 'Color',line_color(5, :), 'LineStyle','--')
end
yline(0, 'LineWidth', 0.8, 'LineStyle', ':')
if strcmp(extra_label, 'q3')
    legend('1st', '2nd', '3rd', 'Orientation', 'horizontal', 'Location', 'north')
elseif strcmp(extra_label, 'q4')
    legend('1st', '2nd', '3rd', '4th', 'Orientation', 'horizontal', ...
        'Location', 'north', 'IconColumnWidth', 28)
elseif strcmp(extra_label, 'q5')
    legend('1st', '2nd', '3rd', '4th', '5th', ...
        'Orientation', 'horizontal', 'Location', 'north', ...
        'IconColumnWidth', 18)
else
    legend('low', 'middle', 'high', 'Orientation', 'horizontal', 'Location', 'north')
end
title('CATE')
xlabel(x_lab)
ylabel('prob. diff.')
xlim(xlim_para)
ylim(ylim_para)
hold off;

for i = 1:n_class
    subplot(plot_rows, 2, i+1)
    plot(x, est_y(:, i), 'LineWidth', 2, 'Color', line_color(i, :))
    hold on;
    CI_plot_X = [x', fliplr(x')];
    CI_plot_Y_95 = [(est_y(:, i) + boot_y_025(:, i))', ...
        fliplr((est_y(:, i) + boot_y_975(:, i))')];
    fill(CI_plot_X, CI_plot_Y_95, [0.8 0.8 0.8], 'FaceAlpha', 0.4, 'EdgeColor', 'none');
    
    CI_plot_Y_90 = [(est_y(:, i) + boot_y_050(:, i))', ...
        fliplr((est_y(:, i) + boot_y_950(:, i))')];
    fill(CI_plot_X, CI_plot_Y_90, [0.65 0.65 0.65], 'FaceAlpha', 0.4, 'EdgeColor', 'none');
    yline(0, 'LineWidth', 0.8, 'LineStyle', ':')
    xlim(xlim_para)
    ylim(ylim_para)
    %title('Low')
    title(title_str{i})
    xlabel(x_lab)
    ylabel('prob. diff.')
    hold off;
end

%%%%%%%%%%
%{
subplot(plot_rows, 2, 2)
plot(x, est_y(:, 1), 'LineWidth', 2, 'Color', [224 143 51]/255)
hold on;
CI_plot_X = [x', fliplr(x')];
CI_plot_Y_95 = [(est_y(:, 1) + boot_y_025(:, 1))', ...
    fliplr((est_y(:, 1) + boot_y_975(:, 1))')];
fill(CI_plot_X, CI_plot_Y_95, [0.8 0.8 0.8], 'FaceAlpha', 0.4, 'EdgeColor', 'none');

CI_plot_Y_90 = [(est_y(:, 1) + boot_y_050(:, 1))', ...
    fliplr((est_y(:, 1) + boot_y_950(:, 1))')];
fill(CI_plot_X, CI_plot_Y_90, [0.65 0.65 0.65], 'FaceAlpha', 0.4, 'EdgeColor', 'none');
yline(0, 'LineWidth', 0.8, 'LineStyle', ':')
xlim(xlim_para)
ylim(ylim_para)
title('Low')
xlabel(x_lab)
ylabel('prob. diff.')
hold off;

subplot(plot_rows, 2, 3)
plot(x, est_y(:, 2), 'LineWidth', 2, 'Color', [246 192 66]/255)
hold on;
CI_plot_X = [x', fliplr(x')];
CI_plot_Y_95 = [(est_y(:, 2) + boot_y_025(:, 2))', ...
    fliplr((est_y(:, 2) + boot_y_975(:, 2))')];
fill(CI_plot_X, CI_plot_Y_95, [0.8 0.8 0.8], 'FaceAlpha', 0.4, 'EdgeColor', 'none');

CI_plot_Y_90 = [(est_y(:, 2) + boot_y_050(:, 2))', ...
    fliplr((est_y(:, 2) + boot_y_950(:, 2))')];
fill(CI_plot_X, CI_plot_Y_90, [0.65 0.65 0.65], 'FaceAlpha', 0.4, 'EdgeColor', 'none');
yline(0, 'LineWidth', 0.8, 'LineStyle', ':')
xlim(xlim_para)
ylim(ylim_para)
title('Middle')
xlabel(x_lab)
ylabel('prob. diff.')
hold off;

subplot(plot_rows, 2, 4)
plot(x, est_y(:, 3),  'LineWidth', 2, 'Color', [170 179 58]/255)
hold on;
CI_plot_X = [x', fliplr(x')];

CI_plot_Y_95 = [(est_y(:, 3) + boot_y_025(:, 3))', ...
    fliplr((est_y(:, 3) + boot_y_975(:, 3))')];
fill(CI_plot_X, CI_plot_Y_95, [0.8 0.8 0.8], 'FaceAlpha', 0.4, 'EdgeColor', 'none');

CI_plot_Y_90 = [(est_y(:, 3) + boot_y_050(:, 3))', ...
    fliplr((est_y(:, 3) + boot_y_950(:, 3))')];
fill(CI_plot_X, CI_plot_Y_90, [0.65 0.65 0.65], 'FaceAlpha', 0.4, 'EdgeColor', 'none');
yline(0, 'LineWidth', 0.8, 'LineStyle', ':')
xlim(xlim_para)
ylim(ylim_para)
title('High')
xlabel(x_lab)
ylabel('prob. diff.')
hold off;
%}

fig = gcf;
fig.Units = 'pixels';
%fig.Position = [406, 229, 674, 420];
fig.Position = [406, 229, 674, 210*plot_rows];
fig.PaperPositionMode = 'auto';

% save plots
print(['./plots/' plot_name '_' extra_label '.eps'], '-depsc2', '-r300')
%print(['./plots/' plot_name '_' extra_label '_10_15.eps'], '-depsc2', '-r300')
end
