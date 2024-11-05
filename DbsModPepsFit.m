 

% This code produces Figure 2C
%NMG 2024 Mayo Clinic

%% load data
load('...path.../LinModelFit.mat');

% Get the field names and initialize variables
fieldNames = fieldnames(LinFit);
slope=nan(length(fieldNames),1);

%% Visualization colors and grouping
ShortLongC = [0.8, 0, 0; 0, 0, .8];         %colors of fit lines, for long vs. short duration DBS
ShortLongCscat = [1, 0.7, .1; .1, .7, 1];   %colors of scatter plot contacts, for long vs. short duration DBS
ShortLongMemb=[1, 1, 2, 1, 1, 1, 2, 1, 2, 1, 2, 2];  %1=long. 2=short. corresponding to pt1L, pt1R, pt2:5, pt6 short, pt6 long, pt7:10

figure
%% Plot scatter points
for i = 1:length(fieldNames)
    % Assuming 'model' is your LinearModel object created using fitlm
    x = LinFit.(fieldNames{i}).Variables.(LinFit.(fieldNames{i}).PredictorNames{1}); % Extract predictor data from model
    y = LinFit.(fieldNames{i}).Variables.(LinFit.(fieldNames{i}).ResponseName); % Extract response data from the same table
    hold on;
    % Plot scatter of original data
    scatter(x, y, 'filled', 'MarkerEdgeColor', 'black', 'MarkerFaceColor', ShortLongCscat(ShortLongMemb(i),:),'MarkerFaceAlpha', 0.5); 
end

%% Plot fit lines and confidence intervals
for i = 1:length(fieldNames)

    % Assuming 'model' is your LinearModel object created using fitlm
    x = LinFit.(fieldNames{i}).Variables.(LinFit.(fieldNames{i}).PredictorNames{1}); % Extract predictor data from model
    y = LinFit.(fieldNames{i}).Variables.(LinFit.(fieldNames{i}).ResponseName); % Extract response data from the same table

    % Generate predictor values for making predictions
    x_fit = linspace(min(x), max(x), 100)'; % Generate a range of x values for fit line

    % Predict the response and confidence intervals
    [y_pred, ci] = predict(LinFit.(fieldNames{i}), x_fit);

    % Plot the regression fit line
    plot(x_fit, y_pred, 'Color', ShortLongC(ShortLongMemb(i),:), 'LineWidth', 4); % 

    % Plot the confidence intervals
    plot(x_fit, ci(:,1), 'Color', ShortLongC(ShortLongMemb(i),:), 'LineStyle', '--', 'LineWidth', 1); % Lower confidence interval, 
    plot(x_fit, ci(:,2), 'Color', ShortLongC(ShortLongMemb(i),:), 'LineStyle', '--', 'LineWidth', 1); % Upper confidence interval, 

    % Enhance the plot
    xlabel(LinFit.(fieldNames{i}).PredictorNames{1});
    ylabel(LinFit.(fieldNames{i}).ResponseName);
    title('Linear Model Fit with Confidence Intervals');
end
hold off;

%% Compare fit slopes short vs long duration DBS
FitSlope = nan(length(fieldNames),1);

for i = 1:length(fieldNames)
    FitSlope(i) = LinFit.(fieldNames{i}).Coefficients{2,1};
end
 
[H,P,ConfI,STATS] = ttest2(FitSlope(ShortLongMemb<1.5), FitSlope(ShortLongMemb>1.5), 'Vartype', 'equal', 'Alpha', 0.05); %

%% Per-subjedt linear regression model metrics in LinFit structure, Supplemental Table 2. e.g.:
% LinFit.pt1L
% 
% ans = % 
% Linear regression model:
%     y ~ 1 + x1
% Estimated Coefficients:
%                    Estimate       SE        tStat       pValue  
%                    ________    ________    _______    __________
% 
%     (Intercept)      1.0669     0.18103     5.8933    6.6523e-08
%     x1             -0.72071    0.082462    -8.7399    1.2974e-13
% 
% Number of observations: 91, Error degrees of freedom: 89
% Root Mean Squared Error: 0.717
% R-squared: 0.462,  Adjusted R-Squared: 0.456
% F-statistic vs. constant model: 76.4, p-value = 1.3e-13