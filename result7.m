% Evaluating the performance of CMIP6 models in simulating Southern Ocean
% biogeochemistry300
% Author: Cheng, M., Ellwood, M., and Maher, N. 
% Last modified: 26/05/2025
% Regression and statistical nalysis (Section 4.1 and 4.2)

all_longname = {'Copernicus','ACCESS-ESM1-5','CanESM5','CESM2','CMCC-ESM2','CNRM-ESM2-1','GFDL-ESM4',...
    'IPSL-CM6A-LR','MPI-ESM-1-2-HAM','MPI-ESM1-2-HR',...
    'MPI-ESM1-2-LR','NorESM2-LM','NorESM2-MM','UKESM1-0-LL'};

mld_djf = [33.23 29.99 38.90 36.71 37.17 25.56 30.82 34.40 20.51 22.33 22.10 42.98 50.52 37.17];
lev_chl = [147.00 90.52 171.26 96.35 152.34 176.51 99.28 150.19 39.38 44.93 43.82 101.27 115.25 126.29];
chl_int = [54.85	26.30	22.06	20.19	26.78	27.39	45.85	30.25	38.84	68.92	34.93	26.44	27.45	36.49];
dcm_fre = [84.71	14.48	84.50	16.52	42.05	70.45	99.63	45.32	4.97	14.50	3.53	0.00	0.00	35.78];

% lev_chl vs mld
p1 = polyfit(mld_djf, lev_chl, 1); 
xfit1 = linspace(min(mld_djf), max(mld_djf), 100);
yfit1 = polyval(p1, xfit1);

% Calculate R²
yfit_all = polyval(p1, mld_djf);
SSres = sum((lev_chl - yfit_all).^2);
SStot = sum((lev_chl - mean(lev_chl)).^2);
R21 = 1 - SSres/SStot;

% Caalculate p calue
X = [ones(length(mld_djf),1) mld_djf'];
[b,~,~,~,stats] = regress(lev_chl', X);
p_val1 = stats(3);

% Plot MLD related trends (Fig.S8)
figs7=figure(1);
set(gcf,'Position',[100 100 1400 600])
subplot(141)
hold on
markers = {'o','s','d','^','v','>','<','p','h','x','+','*','.','diamond'};
colors = lines(length(mld_djf));

for i = 1:length(mld_djf)
    scatter(mld_djf(i), lev_chl(i), 80, 'Marker', markers{mod(i-1,length(markers))+1}, ...
        'MarkerEdgeColor', colors(i,:), 'DisplayName', all_longname{i}, 'LineWidth', 1.5);
end

plot(xfit1, yfit1, 'k--', 'LineWidth', 2, 'DisplayName', 'Regression');

text(18,186.4,'(a)', 'FontWeight', 'bold', 'FontSize', 14);

text(24, 80, sprintf('R^2=%.2f, p=%.3f', R21, p_val1), ...
    'FontSize', 12, 'Color', 'k', 'BackgroundColor', 'none');
xlabel('MLD (m)')
ylabel('Chlorophyll threshold depth (m)')
legend('Location','eastoutside','Position',[0.72 0.1 0.1 0.4])
grid on
box on

% chl_int vs mld
p2 = polyfit(mld_djf, chl_int, 1); 
xfit2 = linspace(min(mld_djf), max(mld_djf), 100);
yfit2 = polyval(p2, xfit2);

yfit_all = polyval(p2, mld_djf);
SSres = sum((chl_int - yfit_all).^2);
SStot = sum((chl_int - mean(chl_int)).^2);
R22 = 1 - SSres/SStot;

X = [ones(length(mld_djf),1) mld_djf'];
[b,~,~,~,stats] = regress(chl_int', X);
p_val2 = stats(3);

subplot(142)
hold on
for i = 1:length(mld_djf)
    scatter(mld_djf(i), chl_int(i), 80, 'Marker', markers{mod(i-1,length(markers))+1}, ...
        'MarkerEdgeColor', colors(i,:), 'LineWidth', 1.5);
end

plot(xfit2, yfit2, 'k--', 'LineWidth', 2);

text(18,72, '(b)', 'FontWeight', 'bold', 'FontSize', 14);

text(28, 42, sprintf('R^2=%.2f, p=%.3f', R22, p_val2), ...
    'FontSize', 12, 'Color', 'k', 'BackgroundColor', 'none');
xlabel('MLD (m)')
ylabel('Integrated chlorophyll in top 100m (mg m^{-2})')
grid on
box on

gaussEqn = 'a*exp(-((x-b)^2)/(2*c^2))';
startPoints = [100, mean(mld_djf), 10]; 
f3 = fit(mld_djf', dcm_fre', gaussEqn, 'Start', startPoints);

xfit3 = linspace(min(mld_djf), max(mld_djf), 200);
yfit3 = f3.a * exp(-((xfit3 - f3.b).^2) / (2 * f3.c^2));

y_pred = f3(mld_djf');
SSres3 = sum((dcm_fre' - y_pred).^2);
SStot3 = sum((dcm_fre' - mean(dcm_fre)).^2);
R23 = 1 - SSres3/SStot3;

subplot(143)
hold on
for i = 1:length(mld_djf)
    scatter(mld_djf(i), dcm_fre(i), 80, 'Marker', markers{mod(i-1,length(markers))+1}, ...
        'MarkerEdgeColor', colors(i,:), 'DisplayName', all_longname{i}, 'LineWidth', 1.5);
end
plot(xfit3, yfit3, 'k--', 'LineWidth', 2, 'DisplayName', 'Gaussian fit');

text(18,104, '(c)', 'FontWeight', 'bold', 'FontSize', 14);

text(37, 60, sprintf('R^2=%.2f', R23), ...
    'FontSize', 12, 'Color', 'k', 'BackgroundColor', 'none');

xlabel('MLD (m)')
ylabel('DCM Frequency (%)')
grid on
box on

print(figs7,'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Figure\figs7','-djpeg','-r600');

% Analyse model components related using in 4.2
% Data
D = [0 0 1 1 1 1 1 0 0 0 0 0 0 1];
T = [0 0 1 0 0 1 0 1 0 0 0 0 0 0];
Reg = [0 1 1 1 1 1 1 0 0 0 0 0 0 0];
Chlv = [0 1 1 1 1 1 1 1 0 0 0 0 0 1];
Fec = [1 NaN 2 2 3 2 3 1 2 2 2 2 2 2];
Siv = [-1 -1 1 1 1 1 1 -1 0 0 0 0 0 1];
KN = [1.00 0.10	0.92 2.25 2.00 2.67	2.00 0.50 0.16 0.16	0.16 0.16 0.16 0.63];
KS = [NaN NaN 0.70 1.00	8.00 2.00 8.00 NaN 1.00	1.00 1.00 5.00 5.00	3.00];
KF = [0.02 NaN 0.18	5.63 6.00 0.37 6.00	0.001 3.60 3.60 3.60 3.60 3.60 0.50];
Chl = [10 11 7 2 9 1 2 5 12 14 13 6 8 4];
NO3 = [6 11 9 10 3 2 1 8 13 14 12 5 7 3];
Si = [NaN NaN 4 8 2 3 1 NaN 10 11 9 5 6 7];
dFe = [5 NaN 11 6 3 8 2 1 11 11 10 7 8 4];
POC = [12 14 8 1 7 3 5 10 4 8 2 10 12 6];
DCMf = [9 1 7 5 2 3 4 10 11 8 12 13 13 6];
model_longname = {'ACCESS-ESM1-5','CanESM5','CESM2','CMCC-ESM2','CNRM-ESM2-1','GFDL-ESM4',...
    'IPSL-CM6A-LR','MIROC-ES2L','MPI-ESM-1-2-HAM','MPI-ESM1-2-HR',...
    'MPI-ESM1-2-LR','NorESM2-LM','NorESM2-MM','UKESM1-0-LL'};
% Chl vs D
Chl_D0 = Chl(D ==0);
Chl_D1 = Chl(D == 1); 

% t-test
[h, p1] = ttest2(Chl_D0, Chl_D1);

% Start plotting (Fig.S9)
fig15=figure(2);
set(gcf,'Position',[100 100 1400 1000])
subplot(3,4,1)
boxplot(Chl,D);
hold on;
for i = 1:length(model_longname)
    scatter(D(i)+1, Chl(i), 40, 'Marker', markers{mod(i-1,length(markers))+1}, ...
        'MarkerEdgeColor', colors(i,:), 'LineWidth', 1.5);
end
xticklabels({'No diatom','Diatom'})
ylabel('Chlorophyll ranking')
ylim([0 14])
set(gca,"YDir",'reverse')
grid on;
text(0.4,-1.5, '(a)', 'FontWeight', 'bold', 'FontSize', 14);
text(2, 12, ['T-test\newlinep = ', num2str(p1, '%.4f')], 'HorizontalAlignment', 'center')

% Chl vs T
Chl_T0 = Chl(T ==0);
Chl_T1 = Chl(T == 1); 

% t-test
[h, p2] = ttest2(Chl_T0, Chl_T1);

subplot(3,4,2)
boxplot(Chl,T);
hold on;
for i = 1:length(model_longname)
    scatter(T(i)+1, Chl(i), 40, 'Marker', markers{mod(i-1,length(markers))+1}, ...
        'MarkerEdgeColor', colors(i,:), 'LineWidth', 1.5);
end
xticklabels({'No diazotroph','Diazotroph'})
ylabel('Chlorophyll ranking')
ylim([0 14])
set(gca,"YDir",'reverse')
grid on;
text(0.4,-1.5, '(b)', 'FontWeight', 'bold', 'FontSize', 14);
text(2, 12, ['T-test\newlinep = ', num2str(p2, '%.4f')], 'HorizontalAlignment', 'center')

% DCMf vs Chlv
DCMf_v0 = DCMf(Chlv ==0);
DCMf_v1 = DCMf(Chlv == 1); 

% t-test
[h, p3] = ttest2(DCMf_v0, DCMf_v1);

subplot(3,4,3)
boxplot(DCMf,Chlv);
hold on;
for i = 1:length(model_longname)
    scatter(Chlv(i)+1, DCMf(i), 40, 'Marker', markers{mod(i-1,length(markers))+1}, ...
        'MarkerEdgeColor', colors(i,:), 'LineWidth', 1.5);
end
xticklabels({'Fixed C:Chl','Variable C:Chl'})
ylabel('DCM frequency ranking')
ylim([0 14])
set(gca,"YDir",'reverse')
grid on;
text(0.4,-1.5, '(c)', 'FontWeight', 'bold', 'FontSize', 14);
text(2, 12, ['T-test\newlinep = ', num2str(p3, '%.4f')], 'HorizontalAlignment', 'center')

% POC vs Siv
% ANOVA
[p4,tbl,stats] = anova1(POC,Siv,'off');

figure(2)
subplot(3,4,5)
boxplot(POC,Siv);
hold on;
for i = 1:length(model_longname)
    scatter(Siv(i)+2, POC(i), 40, 'Marker', markers{mod(i-1,length(markers))+1}, ...
        'MarkerEdgeColor', colors(i,:), 'LineWidth', 1.5);
end
ylabel('POC ranking')
ylim([0 14])
set(gca,"YDir",'reverse','XTickLabel',{'No Si','Fixed C:Si','Variable C:Si'},'FontSize',8)
grid on;
text(0.4,-1.5, '(d)', 'FontWeight', 'bold', 'FontSize', 14);
text(3, 12, ['ANOVA\newlinep = ', num2str(p4, '%.4f')], 'HorizontalAlignment', 'center')

% Si vs KSi
valid_idx = ~isnan(KS) & ~isnan(Si); 
KS_clean = KS(valid_idx);
Si_clean = Si(valid_idx);

lr5 = polyfit(KS_clean, Si_clean, 1); 
xfit5 = linspace(min(KS_clean), max(KS_clean), 100);
yfit5 = polyval(lr5, xfit5);

yfit_all = polyval(lr5, KS_clean);
SSres = sum((Si_clean - yfit_all).^2);
SStot = sum((Si_clean - mean(Si_clean)).^2);
R25 = 1 - SSres/SStot;

X = [ones(length(KS_clean),1) KS_clean'];
[b,~,~,~,stats] = regress(Si_clean', X);
p5 = stats(3);

subplot(3,4,6)
hold on;
for i = 1:length(model_longname)
    scatter(KS(i), Si(i), 40, 'Marker', markers{mod(i-1,length(markers))+1}, ...
        'MarkerEdgeColor', colors(i,:), 'LineWidth', 1.5);
end

plot(xfit5, yfit5, 'k--', 'LineWidth', 2);
xlabel('K_{mSi} (mmol m^{-3})')
ylabel('Silicate ranking')
ylim([0 14])
set(gca,"YDir",'reverse')
grid on;
box on;
text(-0.4,-1.5, '(e)', 'FontWeight', 'bold', 'FontSize', 14);
text(5, 12, sprintf('Linear Regression\nR^2=%.2f, p=%.3f', R25, p5), 'HorizontalAlignment', 'center')

% NO3 vs KN
valid_idx = ~isnan(KN) & ~isnan(NO3); 
KN_clean = KN(valid_idx);
NO3_clean = NO3(valid_idx);

lr5 = polyfit(KN_clean, NO3_clean, 1); 
xfit5 = linspace(min(KN_clean), max(KN_clean), 100);
yfit5 = polyval(lr5, xfit5);

yfit_all = polyval(lr5, KN_clean);
SSres = sum((NO3_clean - yfit_all).^2);
SStot = sum((NO3_clean - mean(NO3_clean)).^2);
R25 = 1 - SSres/SStot;

X = [ones(length(KN_clean),1) KN_clean'];
[b,~,~,~,stats] = regress(NO3_clean', X);
p5 = stats(3);

subplot(3,4,7)
hold on;
for i = 1:length(model_longname)
    scatter(KN(i), NO3(i), 40, 'Marker', markers{mod(i-1,length(markers))+1}, ...
        'MarkerEdgeColor', colors(i,:), 'DisplayName', model_longname{i},'LineWidth', 1.5);
end
plot(xfit5, yfit5, 'k--', 'DisplayName','Regression','LineWidth', 2);
xlabel('K_{mNO3} (mmol m^{-3})')
ylabel('Nitrate ranking')
ylim([0 14])
set(gca,"YDir",'reverse')
grid on;
box on;
text(-0.15,-1.5, '(f)', 'FontWeight', 'bold', 'FontSize', 14);
text(2, 12, sprintf('Linear Regression\nR^2=%.2f, p=%.3f', R25, p5), 'HorizontalAlignment', 'center')
legend('Location','eastoutside','Position',[0.72 0.1 0.1 0.25])

% dFe vs KF
valid_idx = ~isnan(KF) & ~isnan(dFe); 
KF_clean = KF(valid_idx);
dFe_clean = dFe(valid_idx);

lr5 = polyfit(KF_clean, dFe_clean, 1); 
xfit5 = linspace(min(KF_clean), max(KF_clean), 100);
yfit5 = polyval(lr5, xfit5);

yfit_all = polyval(lr5, KF_clean);
SSres = sum((dFe_clean - yfit_all).^2);
SStot = sum((dFe_clean - mean(dFe_clean)).^2);
R25 = 1 - SSres/SStot;

X = [ones(length(KF_clean),1) KF_clean'];
[b,~,~,~,stats] = regress(dFe_clean', X);
p5 = stats(3);

subplot(3,4,9)
hold on;
for i = 1:length(model_longname)
    scatter(KF(i), dFe(i), 40, 'Marker', markers{mod(i-1,length(markers))+1}, ...
        'MarkerEdgeColor', colors(i,:), 'DisplayName', model_longname{i},'LineWidth', 1.5);
end
plot(xfit5, yfit5, 'k--', 'DisplayName','Regression','LineWidth', 2);
xlabel('K_{mFe} (nmol m^{-3})')
ylabel('Dissolved iron ranking')
ylim([0 14])
set(gca,"YDir",'reverse')
grid on;
box on;
text(-0.15,-1.5, '(g)', 'FontWeight', 'bold', 'FontSize', 14);
text(2, 12, sprintf('Linear Regression\nR^2=%.2f, p=%.3f', R25, p5), 'HorizontalAlignment', 'center')
legend('Location','eastoutside','Position',[0.72 0.1 0.1 0.25])

% dFe vs Fec
% ANOVA
[p8,tbl,stats] = anova1(dFe,Fec,'off');

subplot(3,4,10)
boxplot(dFe,Fec);
hold on;
for i = 1:length(model_longname)
    scatter(Fec(i), dFe(i), 40, 'Marker', markers{mod(i-1,length(markers))+1}, ...
        'MarkerEdgeColor', colors(i,:), 'DisplayName', model_longname{i}, 'LineWidth', 1.5);
end
ylabel('Dissolved iron ranking')
ylim([0 14])
set(gca,"YDir",'reverse','XTickLabel',{'Simple','Ligand','Complex'})
grid on;
text(0.4,-1.5, '(h)', 'FontWeight', 'bold', 'FontSize', 14);
text(3, 12, ['ANOVA\newlinep = ', num2str(p8, '%.4f')], 'HorizontalAlignment', 'center')

% DCMf vs Reg
DCMf_r0 = DCMf(Reg ==0);
DCMf_r1 = DCMf(Reg == 1); 

% t-test
[h, p9] = ttest2(DCMf_r0, DCMf_r1);

subplot(3,4,11)
boxplot(DCMf,Reg);
hold on;
for i = 1:length(model_longname)
    scatter(Reg(i)+1, DCMf(i), 40, 'Marker', markers{mod(i-1,length(markers))+1}, ...
        'MarkerEdgeColor', colors(i,:), 'LineWidth', 1.5);
end
xticklabels({'No ammonium','Uptake ammonium'})
ylabel('DCM frequency ranking')
ylim([0 14])
set(gca,"YDir",'reverse')
grid on;
text(0.4,-1.5, '(i)', 'FontWeight', 'bold', 'FontSize', 14);
text(2, 12, ['T-test\newlinep = ', num2str(p9, '%.4f')], 'HorizontalAlignment', 'center')

print(fig15,'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Figure\fig15','-djpeg','-r600');
