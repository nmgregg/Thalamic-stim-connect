% 
% This code produces Figure 4
% 
% NMG 2024 Mayo Clinic

clear all

%% load data

load('...path.../PepsIEDsData.mat')

%% resize to 1d arrays
xr=group.RMSCohenD'; %reshape 1D  (subj x channels, unwrapped into 1d)
xr=xr(:);
yrtempB=group.SpikesBase'; %same
yrtempB=yrtempB(:);
yrtempOn=group.SpikesStimOn'; %same
yrtempOn=yrtempOn(:);
yr=yrtempOn./yrtempB;

%% remove nan vals
validIndices = ~isnan(xr) & ~isnan(yr);
xr = xr(validIndices);
yr = yr(validIndices);

%% Initialize and Calculate Average Values for Seizure Network (SN)
group.RMSCohenDSzNetOnly=nan(size(group.RMSCohenD)); %initialize variable; will only keep seizure network PEPs
group.RMSCohenDSzNetAvg=nan(size(group.RMSCohenD));  %initialize variable; replace each SN contact value with average SN PEP
for i=1:length(group.subj)
    group.RMSCohenDSzNetOnly(i, isfinite(group.SpikesBase(i,:)))=group.RMSCohenD(i, isfinite(group.SpikesBase(i,:)));
    group.RMSCohenDSzNetAvg(i,:)=nanmean(group.RMSCohenDSzNetOnly(i,:),2);
end

groupRMSCohenDSzNetAvg=reshape(group.RMSCohenDSzNetAvg', numel(group.RMSCohenDSzNetAvg), 1);
xr_AvgSNvals=groupRMSCohenDSzNetAvg(validIndices); %keep only valid SN values. Use yr for spike rate values, these are not averaged across SN contacts for this variable.


%% Define and fit exponential model
expDecayModel = fittype('a*exp(-b*x)+c', 'independent', 'x', 'coefficients', {'a', 'b', 'c'});
[fitResultr, gofr] = fit(xr, yr, expDecayModel);

%% Display fitting results
disp('Fitted Model:');
disp(fitResultr);
disp('Goodness of Fit:');
disp(gofr);

%% Model and Confidence Intervals
coeffsr = coeffvalues(fitResultr);
confIntsr = confint(fitResultr);
numCoeffsr = numel(coeffsr);

fprintf('Coefficient Estimates and 95%% Confidence Intervals:\n');
for i = 1:numCoeffsr
    fprintf('Coefficient %d: %f (%f, %f)\n', i, coeffsr(i), confIntsr(1, i), confIntsr(2, i));
end

%% Plot fit with data

figure(101);

subplot(1,2,1)
x_fitr = linspace(min(xr), max(xr), 100);
y_fitr = feval(fitResultr, x_fitr);
cir = predint(fitResultr, x_fitr, 0.95, 'functional', 'off');

plot(xr, yr, 'bo');
hold on;
plot(x_fitr, y_fitr, 'r-');
fill([x_fitr, fliplr(x_fitr)], [cir(:,1); flipud(cir(:,2))], 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
xlabel('x');
ylabel('y');
title('Fit with 95% Confidence Intervals');
legend('Data', 'Fit', '95% Confidence Interval');
hold off;

%% Stat sig of fit for each coefficient

nr = length(yr); % number of data points
pr = numCoeffsr; % number of coefficients

% Degrees of freedom
dofr = nr - pr;

% Standard errors of the coefficients
ser = (confIntsr(2, :) - confIntsr(1, :)) / (2 * 1.96);

% t-statistics for the coefficients
tstatsr = coeffsr ./ ser;

% p-values for the coefficients
pvalsr = 2 * (1 - tcdf(abs(tstatsr), dofr));

fprintf('P-values for the coefficients:\n');
for i = 1:numCoeffsr
    fprintf('Coefficient %d: p-value = %f\n', i, pvalsr(i)); %
end

%% Analysis for average SN vals
spikemod=group.SpikesStimOn./group.SpikesBase; %relative spike rate with stim vs. baseline

for i=1:length(group.subj)
    xav_r(i)=nanmean(group.RMSCohenD(i,isfinite(spikemod(i,:))),2);
    yav_r(i)=nanmean(spikemod(i,isfinite(spikemod(i,:))), 2);

end
xav_r=xav_r'; yav_r=yav_r';
[fitResultav_r, gofav_r, outputav_r] = fit(xav_r, yav_r, expDecayModel);


%% display

disp('Fitted Model:');
disp(fitResultav_r);
disp('Goodness of Fit:');
disp(gofav_r);

%% Model and Confidence Intervals

coeffsav_r = coeffvalues(fitResultav_r);
confIntsav_r = confint(fitResultav_r);
numCoeffsav_r = numel(coeffsav_r);

fprintf('Coefficient Estimates and 95%% Confidence Intervals:\n');
for i = 1:numCoeffsav_r
    fprintf('Coefficient %d: %f (%f, %f)\n', i, coeffsav_r(i), confIntsav_r(1, i), confIntsav_r(2, i));
end


%% Stat sig of fit

nav_r = length(yav_r); % number of data points
pav_r = numCoeffsav_r; % number of coefficients

% Degrees of freedom
dofav_r = nav_r - pav_r;

% Standard errors of the coefficients
seav_r = (confIntsav_r(2, :) - confIntsav_r(1, :)) / (2 * 1.96);

% t-statistics for the coefficients
tstatsav_r = coeffsav_r ./ seav_r;

% p-values for the coefficients
pvalsav_r = 2 * (1 - tcdf(abs(tstatsav_r), dofav_r));

fprintf('P-values for the coefficients:\n');
for i = 1:numCoeffsav_r
    fprintf('SN avg Coefficient %d: p-value = %f\n', i, pvalsr(i));
end

%% define colors for 11 thalami 
colors = [
    0.0000, 0.4470, 0.7410;  % blue
    0.8500, 0.3250, 0.0980;  % orange
    0.9290, 0.6940, 0.1250;  % yellow
    0.4940, 0.1840, 0.5560;  % purple
    0.4660, 0.6740, 0.1880;  % green
    0.3010, 0.7450, 0.9330;  % light blue
    0.6350, 0.0780, 0.1840;  % dark red
    0.0000, 0.5000, 0.0000;  % green
    1.0000, 0.0000, 0.0000;  % red
    0.0000, 1.0000, 1.0000;  % cyan
    1.0000, 0.0000, 1.0000;  % magenta
];

%%
figure(101)
subplot(1,2,2)
x_fitav_r = linspace(min(xav_r), max(xav_r), 100);
y_fitav_r = feval(fitResultav_r, x_fitav_r);
ciav_r = predint(fitResultav_r, x_fitav_r, 0.95, 'functional', 'off');

plot(xav_r, yav_r, 'bo');
hold on;
plot(x_fitav_r, y_fitav_r, 'r-');
fill([x_fitav_r, fliplr(x_fitav_r)], [ciav_r(:,1); flipud(ciav_r(:,2))], 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
xlabel('x');
ylabel('y');
title('Fit with 95% Confidence Intervals');
legend('Data', 'Fit', '95% Confidence Interval');
hold off;
xlabel(['r^2=',num2str(gofav_r.rsquare)])
ylabel(['fit model: fitResult(x) = a*exp(-b*x)+c'])

%%

for i=1:11
    subplot(1,2,1); hold on; 
    plot(group.RMSCohenD(i,:), spikemod(i,:),'.','MarkerEdgeColor', colors(i,:), 'MarkerSize',20); 
 end
xlabel(['r^2=',num2str(gofr.rsquare)])
legend(['data','fit','95%CI',group.subj]); % or use 'Location', 'northeast' for a fixed position

for i=1:11
    subplot(1,2,2); hold on; 
    plot(xav_r(i), yav_r(i),'.','MarkerEdgeColor', colors(i,:), 'MarkerSize',40);
end
xlabel(['r^2=',num2str(gofav_r.rsquare)])
legend(['data','fit','95%CI',group.subj]); % or use 'Location', 'northeast' for a fixed position

%% Binarize and boxplot, for seizure network chans with baseline PEP effect size >1 vs. < 1
groupBin=nan(length(xr),2);
groupBin(xr<1,1)=yr(xr<1);
groupBin(xr>1,2)=yr(xr>1);
[P,ANOVATAB,STATS] = anova1(groupBin);
hold on;
markerSize=40;
scatter(ones(length(groupBin(:,1)),1), groupBin(:,1), markerSize, 'filled', 'jitter', 'on', 'jitterAmount', 0.15);
scatter(2*ones(length(groupBin(:,1)),1), groupBin(:,2), markerSize, 'filled', 'jitter', 'on', 'jitterAmount', 0.15);

%% Binarize and boxplot, for seizure network average PEP effect size >0.6 vs <0.6 
groupAvgBin=nan(length(xr_AvgSNvals),2);
groupAvgBin(xr_AvgSNvals<0.6,1)=yr(xr_AvgSNvals<0.6);
groupAvgBin(xr_AvgSNvals>0.6,2)=yr(xr_AvgSNvals>0.6);
[Pavg,ANOVATABavg,STATSavg] = anova1(groupAvgBin);
hold on;
markerSize=40;
scatter(ones(length(groupAvgBin(:,1)),1), groupAvgBin(:,1), markerSize, 'filled', 'jitter', 'on', 'jitterAmount', 0.15);
scatter(2*ones(length(groupAvgBin(:,1)),1), groupAvgBin(:,2), markerSize, 'filled', 'jitter', 'on', 'jitterAmount', 0.15);
