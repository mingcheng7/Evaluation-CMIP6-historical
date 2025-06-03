% Evaluating the performance of CMIP6 models in simulating Southern Ocean
% biogeochemistry
% Author: Cheng, M., Ellwood, M., and Maher, N. 
% Last modified: 26/05/2025
% Heat map for ranking (Section 3.4)

% Load data
rankfile = 'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Data\Ranking.txt';
rawdata = readtable(rankfile,'Delimiter','\t');
rank1 = table2array(rawdata(:,1:7));
rank2 = table2array(rawdata(:,8:14));
varnames = {'chl','NO3','Si','dFe','POC','DCM','OVR'};
model_longname = {'ACCESS-ESM1-5','CanESM5','CESM2','CMCC-ESM2','CNRM-ESM2-1','GFDL-ESM4',...
    'IPSL-CM6A-LR','MIROC-ES2L','MPI-ESM-1-2-HAM','MPI-ESM1-2-HR',...
    'MPI-ESM1-2-LR','NorESM2-LM','NorESM2-MM','UKESM1-0-LL'};

% Plot heat map (Fig.13)
fig13=figure(1);
set(gcf,'Position',[100 100 760 1000])
imagesc(rank1,'AlphaData',~isnan(rank1));
xticklabels(varnames)
yticks(1:14)
yticklabels(model_longname)
colormap(flipud(nclCM(143,14)))
set(gca,'Color',[.5 .5 .5]);
clim([0.5 14.5])
set(gca,'XAxisLocation','top')
for i = 1:size(rank1, 1)
    for j = 1:size(rank1, 2)
        value = rank1(i, j);
        if ~isnan(value)  
            text(j, i, sprintf('%.0f', value), ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle', ...
                'Color', 'k', ...
                'FontSize', 12);
        end
    end
end
ax = gca;
ax.XAxis.FontSize = 14; 
ax.YAxis.FontSize = 14;

print(fig13,'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Figure\fig13','-djpeg','-r300');

