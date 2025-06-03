% Evaluating the performance of CMIP6 models in simulating Southern Ocean
% biogeochemistry
% Author: Cheng, M., Ellwood, M., and Maher, N. 
% Last modified: 26/05/2025
% Dissolved iron related analysis (Section 3.1)
clear;

% Load data
dfefile_obs = 'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\In-situ\tagliabue_fe_database_jun2015_public.xlsx';
dfefile_acc = 'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Regrided_data\CMIP\dfe_ACCESS-ESM1-5_final.nc';
dfefile_ces = 'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Regrided_data\CMIP\dfe_CESM2_final.nc';
dfefile_cmc = 'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Regrided_data\CMIP\dfe_CMCC-ESM2_final.nc';
dfefile_cnr = 'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Regrided_data\CMIP\dfe_CNRM-ESM2-1_final.nc';
dfefile_gfd = 'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Regrided_data\CMIP\dfe_GFDL-ESM4_final.nc';
dfefile_ips = 'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Regrided_data\CMIP\dfe_IPSL-CM6A-LR_final.nc';
dfefile_mir = 'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Regrided_data\CMIP\dfe_MIROC-ES2L_final.nc';
dfefile_mhm = 'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Regrided_data\CMIP\dfe_MPI-ESM-1-2-HAM_final.nc';
dfefile_mhr = 'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Regrided_data\CMIP\dfe_MPI-ESM1-2-HR_final.nc';
dfefile_mlr = 'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Regrided_data\CMIP\dfe_MPI-ESM1-2-LR_final.nc';
dfefile_nlm = 'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Regrided_data\CMIP\dfe_NorESM2-LM_final.nc';
dfefile_nmm = 'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Regrided_data\CMIP\dfe_NorESM2-MM_final.nc';
dfefile_uke = 'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Regrided_data\CMIP\dfe_UKESM1-0-LL_final.nc';

model_abbs = {'acc','ces','cmc','cnr','gfd','ips','mir','mhm','mhr','mlr','nlm','nmm','uke'};
model_longname = {'ACCESS-ESM1-5','CESM2','CMCC-ESM2','CNRM-ESM2-1','GFDL-ESM4',...
    'IPSL-CM6A-LR','MIROC-ES2L','MPI-ESM-1-2-HAM','MPI-ESM1-2-HR',...
    'MPI-ESM1-2-LR','NorESM2-LM','NorESM2-MM','UKESM1-0-LL'};
lev_name = {'st_ocean','lev','lev','lev','lev','olevel','lev','lev','lev','lev','lev','lev','lev'};

dfe_obs_data = readtable(dfefile_obs,'Sheet',2);
dfe_obs_mat = table2array(dfe_obs_data(:,1:6));
dfe_obs_mat((dfe_obs_mat(:,6)<0),6) = NaN;
ind_os = (dfe_obs_mat(:,5)<=30)&(dfe_obs_mat(:,4)<=-30);
mon_obs = dfe_obs_mat(ind_os,1);
yrs_obs = dfe_obs_mat(ind_os,2);
lon_obs = dfe_obs_mat(ind_os,3);
lat_obs = dfe_obs_mat(ind_os,4);
lev_obs = dfe_obs_mat(ind_os,5);
dfe_obs = dfe_obs_mat(ind_os,6);
lon_obsp = lon_obs;
lon_obs(lon_obs<0) = lon_obs(lon_obs<0) + 360;

lon = ncread(dfefile_ces,'lon');
lat = ncread(dfefile_ces,'lat');

lonr = [lon;lon(1)];

[lat2,lon2] = meshgrid(lat,lon);
[lat2r,lon2r] = meshgrid(lat,lonr);

for i = 1:length(model_abbs)
    filename = eval(['dfefile_' model_abbs{i}]);
    dfe_read = ncread(filename,'dfe');
    lev_read = ncread(filename,lev_name{i});
    eval(['dfe_' model_abbs{i} '_read=dfe_read;']);
    eval(['lev_' model_abbs{i} '=lev_read;']);
end

% Optimise some depth mannually
lev = [5;15;25;35;45;55;65;75;85;100;125;150;175;200];
lev_mhm(1) = 5;
lev_mhr(1) = 5;
lev_mlr(1) = 5;

% Regrid dfe data
[~,lon_idx] = histc(lon_obs,lon);
[~,lat_idx] = histc(lat_obs,lat);

valid_idx = lon_idx > 0 & lon_idx < length(lon) & ...
            lat_idx > 0 & lat_idx < length(lat);
lon_idx = lon_idx(valid_idx);
lat_idx = lat_idx(valid_idx);
dfe_obs2 = dfe_obs(valid_idx);

grid_idx = sub2ind([360,60],lon_idx,lat_idx);

sum_dfe = accumarray(grid_idx,dfe_obs2,[360*60,1],@sum,NaN);
count_dfe = accumarray(grid_idx,1,[360*60,1],@sum,NaN);
avg_dfe = sum_dfe./count_dfe;

dfe_obs_grid = reshape(avg_dfe,[360,60]);

% Interplate CMIP6 data to (virtual) real data grid size
for i = 1:length(model_abbs)
    dfe_regrid = NaN(360,60,length(lev),180);
    dfe_read = eval(['dfe_' model_abbs{i} '_read']);
    lev_read = eval(['lev_' model_abbs{i}]);
    for ji = 1:size(dfe_regrid,1)
        for jj = 1:size(dfe_regrid,2)
            for jk = 1:size(dfe_regrid,4)
                dfe_regrid(ji,jj,:,jk) = interp1(lev_read,squeeze(dfe_read(ji,jj,:,jk)),lev)*1e6;
            end
        end
    end
    dfe361 = dfe_regrid([1:360,1],:,:,:);
    eval(['dfe_' model_abbs{i} '=dfe361;']);
end

dfe_acc = dfe_acc*1e-6;

% Fill blank for CMCC-ESM2, CanESM5, IPSL-CM6A-LR, MPI-ESM-1-2-HAM
dfe_cmc_regrid = dfe_cmc;
dfe_cnr_regrid = dfe_cnr;
dfe_ips_regrid = dfe_ips;
dfe_mhm_regrid = dfe_mhm;
dfe_mlr_regrid = dfe_mlr;

dfe_cmc(73,:,:,:) = mean(dfe_cmc_regrid([72,75],:,:,:));
dfe_cmc(74,:,:,:) = mean(dfe_cmc_regrid([72,75],:,:,:));
dfe_cnr(73,:,:,:) = mean(dfe_cnr_regrid([72,75],:,:,:));
dfe_cnr(74,:,:,:) = mean(dfe_cnr_regrid([72,75],:,:,:));
dfe_ips(73,:,:,:) = mean(dfe_ips_regrid([72,75],:,:,:));
dfe_ips(74,:,:,:) = mean(dfe_ips_regrid([72,75],:,:,:));

lon_fill_mhm = [147;148;148;149;149;150;150;151;151;152;152;153;153;153;...
    154;154;154;155;155;155;155;156;156;156;156;157;157;157;157;157];
lat_fill_mhm = [24;25;26;27;28;29;30;31;32;34;35;37;38;39;41;42;43;45;46;...
    47;48;50;51;52;53;56;57;58;59;60];

for f = 1:length(lon_fill_mhm)
    for ji = 1:length(lon_fill_mhm)
        for jj = 1:size(dfe_mhm,3)
            for jk = 1:size(dfe_mhm,4)
                dfe_mhm(lon_fill_mhm(ji),lat_fill_mhm(ji),jj,jk) = ...
                    mean(dfe_mhm_regrid([lon_fill_mhm(ji)-1,lon_fill_mhm(ji)+1],lat_fill_mhm(ji),jj,jk));
            end
        end
    end
end

for f = 1:length(lon_fill_mhm)
    for ji = 1:length(lon_fill_mhm)
        for jj = 1:size(dfe_mlr,3)
            for jk = 1:size(dfe_mlr,4)
                dfe_mlr(lon_fill_mhm(ji),lat_fill_mhm(ji),jj,jk) = ...
                    mean(dfe_mlr_regrid([lon_fill_mhm(ji)-1,lon_fill_mhm(ji)+1],lat_fill_mhm(ji),jj,jk));
            end
        end
    end
end

% Calculate mean dfe in austral summer (DJJ) from 2000 to 2014
idx = ([1,2,12].*ones(15,1)+(0:12:14*12)')';
idx1 = idx(:);

for i = 1:length(model_abbs)
    dfe_regrid = eval(['dfe_' model_abbs{i}]);
    dfe_djj_mean = mean(dfe_regrid(:,:,:,idx1),4);
    eval(['dfe_' model_abbs{i} '_djj_mean=dfe_djj_mean;']);
end

% Interplate surface dfe CMIP6 data to real grid
for i = 1:length(model_abbs)
    dfe_mod1 = NaN(size(dfe_obs));
    dfe_djj_mean = eval(['dfe_' model_abbs{i} '_djj_mean']);
    lev_read = eval(['lev_' model_abbs{i}]);
    for j = 1:length(dfe_mod1)
        dfe_mod1(j) = interp2(lat2,lon2,dfe_djj_mean(1:360,:,1),lat_obs(j),lon_obs(j));
    end
    dfe_diff = dfe_mod1 - dfe_obs;
    eval(['dfe_' model_abbs{i} '_os=dfe_mod1;']);
    eval(['dfe_' model_abbs{i} '_diff_os=dfe_diff;']);
end

% Read front data
frontname = {'stf','saf','pf'};
for i = 1:length(frontname)
    txtfile = ['D:\OneDrive - Australian National University\PhD\Project1\MidTerm\Data\' frontname{i} '.txt'];
    fid = fopen(txtfile,'r');
    header = fgetl(fid);
    data = textscan(fid,'%f %f');
    fclose(fid);
    line=[];
    line(:,1) = data{1};
    line(:,2) = data{2};
    eval(['line_' frontname{i} '=line;']);
end

% Define subplot number
splist = [2;3;4;5;6;8;9;10;11;12;14;15;16;17];

% Plot surface dissolved iron (Fig.7)
fig7=figure(1);
set(gcf,'Position',[100 100 1200 600])
sp0=subplot(3,6,1);
pos = get(gca, 'Position'); 
set(gca, 'Position', [pos(1) pos(2) pos(3)*1.15 pos(4)*1.15]); 
m_proj('stereographic','lat',-90,'long',0,'radius',60);
m_scatter(lon_obs,lat_obs,20,dfe_obs,'o','filled');
shading flat;
colormap(sp0,flipud(m_colmap('blues')));
m_line(line_stf(:,1),line_stf(:,2),'color','k','LineStyle','--','linewi',.5)
m_line(line_saf(:,1),line_saf(:,2),'color','k','LineStyle','--','linewi',.5)
m_line(line_pf(:,1),line_pf(:,2),'color','k','LineStyle','--','linewi',.5)
m_coast('patch','w');
m_coast('linewidth',1,'color','k');
m_grid('xtick',6,'xticklabel',[],'tickdir','out','ytick',[],'linest','none','xaxisloc','top','yaxislocation','middle','backcolor',[0.7 0.7 0.7]);
title('GEOTRACERS','FontSize',12)
clim([0 1])
c=colorbar;
c.Ticks = 0:0.25:1;
c.Label.String = 'dFe (\mumol/m^{3})';
c.Label.FontSize = 11;
set(c,'Location','southoutside','Position',[0.14, 0.65, 0.1, 0.02])

dfe_can_diff_os = NaN(size(dfe_obs));
plot_abbs = {'acc','can','ces','cmc','cnr','gfd','ips','mir','mhm','mhr','mlr','nlm','nmm','uke'};
plot_longname = {'ACCESS-ESM1-5','CanESM5^{*}','CESM2','CMCC-ESM2','CNRM-ESM2-1','GFDL-ESM4',...
    'IPSL-CM6A-LR','MIROC-ES2L','MPI-ESM-1-2-HAM','MPI-ESM1-2-HR',...
    'MPI-ESM1-2-LR','NorESM2-LM','NorESM2-MM','UKESM1-0-LL'};

for s = 1:length(splist)
    plotdata = eval(['dfe_' plot_abbs{s} '_diff_os']);
    sp=subplot(3,6,splist(s));
    pos = get(gca, 'Position'); 
    set(gca, 'Position', [pos(1) pos(2) pos(3)*1.15 pos(4)*1.15]); 
    m_proj('stereographic','lat',-90,'long',0,'radius',60);
    m_scatter(lon_obs,lat_obs,20,plotdata,'o','filled');
    shading flat;
    colormap(sp,b2r(-1,1));
    m_line(line_stf(:,1),line_stf(:,2),'color','k','LineStyle','--','linewi',.5)
    m_line(line_saf(:,1),line_saf(:,2),'color','k','LineStyle','--','linewi',.5)
    m_line(line_pf(:,1),line_pf(:,2),'color','k','LineStyle','--','linewi',.5)
    m_coast('patch','w');
    m_coast('linewidth',1,'color','k');
    m_grid('xtick',6,'xticklabel',[],'tickdir','out','ytick',[],'linest','none','xaxisloc','top','yaxislocation','middle','backcolor',[0.7 0.7 0.7]);
    title(plot_longname{s},'FontSize',12)
end
c1=colorbar;
set(c1,'Location','eastoutside','Position',[0.82, 0.12, 0.01, 0.25])
c1.Label.String = 'dFe (\mumol/m^{3})';
c1.Label.FontSize = 11;

print(fig7,'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Figure\fig7','-djpeg','-r300');

% Calculate mean bias error for 4 zones and total SO
so_boundary1 = [-180,-30;180,-30];
so_boundary2 = [-180,-90;180,-90];
poly_stz = [so_boundary1;flipud(line_stf)];
poly_saz = [line_stf;flipud(line_saf)];
poly_pfz = [line_saf;flipud(line_pf)];
poly_az = [line_pf;flipud(so_boundary2)];

idx_stz = inpolygon(lon_obsp,lat_obs,poly_stz(:,1),poly_stz(:,2));
idx_saz = inpolygon(lon_obsp,lat_obs,poly_saz(:,1),poly_saz(:,2));
idx_pfz = inpolygon(lon_obsp,lat_obs,poly_pfz(:,1),poly_pfz(:,2));
idx_az = inpolygon(lon_obsp,lat_obs,poly_az(:,1),poly_az(:,2));

mbe_dfe_os = nan(length(plot_abbs),5);
for i = 1:length(plot_abbs)
    dfe_diff_os = eval(['dfe_' plot_abbs{i} '_diff_os(:,:,1)']);
    mbe_dfe_os(i,1) = sum(dfe_diff_os(~isnan(dfe_diff_os)).*cosd(lat_obs(~isnan(dfe_diff_os))))...
        ./sum((cosd(lat_obs(~isnan(dfe_diff_os)))));
    mbe_dfe_os(i,2) = sum(dfe_diff_os(~isnan(dfe_diff_os)&idx_stz).*cosd(lat_obs(~isnan(dfe_diff_os)&idx_stz)))...
        ./sum(cosd(lat_obs(~isnan(dfe_diff_os)&idx_stz)));
    mbe_dfe_os(i,3) = sum(dfe_diff_os(~isnan(dfe_diff_os)&idx_saz).*cosd(lat_obs(~isnan(dfe_diff_os)&idx_saz)))...
        ./sum(cosd(lat_obs(~isnan(dfe_diff_os)&idx_saz)));
    mbe_dfe_os(i,4) = sum(dfe_diff_os(~isnan(dfe_diff_os)&idx_pfz).*cosd(lat_obs(~isnan(dfe_diff_os)&idx_pfz)))...
        ./sum(cosd(lat_obs(~isnan(dfe_diff_os)&idx_pfz)));
    mbe_dfe_os(i,5) = sum(dfe_diff_os(~isnan(dfe_diff_os)&idx_az).*cosd(lat_obs(~isnan(dfe_diff_os)&idx_az)))...
        ./sum(cosd(lat_obs(~isnan(dfe_diff_os)&idx_az)));
end

% Plot MBE (Fig.8)
col = [0.4 0 0.6; 1 0.7 0.7; 0 0.8 0; 1 0.8 0; 0.4 0.7 1];
fig8=figure(2);
set(gcf,'Position',[100,100,1200,500])
b=bar(mbe_dfe_os);
for c = 1:length(b)
    b(c).FaceColor = col(c,:);
end
ylabel('Mean Bias Error (\mumol/m^{3})')
set(gca,'XTickLabel',plot_longname,'FontSize',12)
legend({'SO','STZ','SAZ','PFZ','AZ'},'Location','southeast')

print(fig8,'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Figure\fig8','-djpeg','-r300');

% Taylor Diagram
% Calculation of different parameters
tdlist = {'obs','acc','can','ces','cmc','cnr','gfd','ips','mir','mhm','mhr','mlr','nlm','nmm','uke'};
dfe_obs_os = dfe_obs;
dfe_can_os = NaN(size(dfe_obs_os));
for i = 1:length(tdlist)
    dfe_so = eval(['dfe_' tdlist{i} '_os']);
    dfe_stz = dfe_so(idx_stz);
    dfe_saz = dfe_so(idx_saz);
    dfe_pfz = dfe_so(idx_pfz);
    dfe_az = dfe_so(idx_az);
    eval(['dfe_so_' tdlist{i} '=dfe_so(:);']);
    eval(['dfe_stz_' tdlist{i} '=dfe_stz;']);
    eval(['dfe_saz_' tdlist{i} '=dfe_saz;']);
    eval(['dfe_pfz_' tdlist{i} '=dfe_pfz;']);
    eval(['dfe_az_' tdlist{i} '=dfe_az;']);
end

rglist = {'so','stz','saz','pfz','az'};

valid_idx_so = ~isnan(dfe_so_obs) & ~isnan(dfe_so_ces) & ...
    ~isnan(dfe_so_cmc) & ~isnan(dfe_so_gfd) & ~isnan(dfe_so_ips) & ...
    ~isnan(dfe_so_mhm) & ~isnan(dfe_so_mhr) & ~isnan(dfe_so_mlr) & ...
    ~isnan(dfe_so_nlm) & ~isnan(dfe_so_nmm);

valid_idx_stz = ~isnan(dfe_stz_obs) & ~isnan(dfe_stz_ces) & ...
    ~isnan(dfe_stz_cmc) & ~isnan(dfe_stz_gfd) & ~isnan(dfe_stz_ips) & ...
    ~isnan(dfe_stz_mhm) & ~isnan(dfe_stz_mhr) & ~isnan(dfe_stz_mlr) & ...
    ~isnan(dfe_stz_nlm) & ~isnan(dfe_stz_nmm);

valid_idx_saz = ~isnan(dfe_saz_obs) & ~isnan(dfe_saz_ces) & ...
    ~isnan(dfe_saz_cmc) & ~isnan(dfe_saz_gfd) & ~isnan(dfe_saz_ips) & ...
    ~isnan(dfe_saz_mhm) & ~isnan(dfe_saz_mhr) & ~isnan(dfe_saz_mlr) & ...
    ~isnan(dfe_saz_nlm) & ~isnan(dfe_saz_nmm);

valid_idx_pfz = ~isnan(dfe_pfz_obs) & ~isnan(dfe_pfz_ces) & ...
    ~isnan(dfe_pfz_cmc) & ~isnan(dfe_pfz_gfd) & ~isnan(dfe_pfz_ips) & ...
    ~isnan(dfe_pfz_mhm) & ~isnan(dfe_pfz_mhr) & ~isnan(dfe_pfz_mlr) & ...
    ~isnan(dfe_pfz_nlm) & ~isnan(dfe_pfz_nmm);

valid_idx_az = ~isnan(dfe_az_obs) & ~isnan(dfe_az_ces) & ...
    ~isnan(dfe_az_cmc) & ~isnan(dfe_az_gfd) & ~isnan(dfe_az_ips) & ...
    ~isnan(dfe_az_mhm) & ~isnan(dfe_az_mhr) & ~isnan(dfe_az_mlr) & ...
    ~isnan(dfe_az_nlm) & ~isnan(dfe_az_nmm);

for j = 1:length(rglist)
    for i = 1:length(tdlist)
        eval(['dfe_' rglist{j} '_' tdlist{i} '_nonan=dfe_' rglist{j} '_' tdlist{i} '(valid_idx_' rglist{j} ');']);
    end
end

stats = [];
for j = 1:length(rglist)
    for i = 1:length(plot_abbs)
        S = SStats(eval(['dfe_' rglist{j} '_obs_nonan']),eval(['dfe_' rglist{j} '_' plot_abbs{i} '_nonan']));
        stats(j,:,i+1) = S;
    end
    stats(j,:,1) = SStats(eval(['dfe_' rglist{j} '_obs_nonan']),eval(['dfe_' rglist{j} '_obs_nonan']));

    stdev = stats(j,2,1); 

    stats(j,2,:) = stats(j,2,:)/stdev; 
    stats(j,3,:) = stats(j,3,:)/stdev;
end

stats(:,1,:) = [];

% Taylor Diagram for surface dissolved iron (Fig.S4)
td_longname = {' ','ACCESS-ESM1-5','CanESM5^{*}','CESM2','CMCC-ESM2','CNRM-ESM2-1','GFDL-ESM4',...
    'IPSL-CM6A-LR','MIROC-ES2L','MPI-ESM-1-2-HAM','MPI-ESM1-2-HR',...
    'MPI-ESM1-2-LR','NorESM2-LM','NorESM2-MM','UKESM1-0-LL'};

figs4=figure(3);
set(gcf,'Position',[100 100 1500 1000])
subplot(2,3,1)
[hp, ht, axl] = taylor_diagram(squeeze(stats(1,1,:)),squeeze(stats(1,2,:)),squeeze(stats(1,3,:)), ...
    'markerLabel',td_longname,'markerLegend', 'on', ...
    'styleSTD', '-', 'colOBS','k', 'markerObs','o', ...
    'markerSize',6, 'tickRMS',0:0.5:1, ...
    'tickRMSangle', 115, 'showlabelsRMS', 'on', ...
    'titleRMS','on', 'titleOBS','Obs');
t1=title('Surface dissolved Fe in SO','FontSize',14);
t1.Position(2) = t1.Position(2) + 0.3;

subplot(2,3,2)
[hp, ht, axl] = taylor_diagram(squeeze(stats(2,1,:)),squeeze(stats(2,2,:)),squeeze(stats(2,3,:)), ...
    'markerLabel',td_longname, 'markerLegend', 'on', ...
    'styleSTD', '-', 'colOBS','k', 'markerObs','o', ...
    'markerSize',6, 'tickRMS',0:0.5:1, ...
    'tickRMSangle', 115, 'showlabelsRMS', 'on', ...
    'titleRMS','on', 'titleOBS','Obs');
t2=title('Surface dissolved Fe in STZ','FontSize',14);
t2.Position(2) = t2.Position(2) + 0.1;

subplot(2,3,3)
[hp, ht, axl] = taylor_diagram(squeeze(stats(3,1,:)),squeeze(stats(3,2,:)),squeeze(stats(3,3,:)), ...
    'markerLabel',td_longname, 'markerLegend', 'on', ...
    'styleSTD', '-', 'colOBS','k', 'markerObs','o', ...
    'markerSize',6, 'tickRMS',0:0.5:1, ...
    'tickRMSangle', 115, 'showlabelsRMS', 'on', ...
    'titleRMS','on', 'titleOBS','Obs');
t3=title('Surface dissolved Fe in SAZ','FontSize',14);
t3.Position(2) = t3.Position(2) + 0.3;

subplot(2,3,4)
[hp, ht, axl] = taylor_diagram(squeeze(stats(4,1,:)),squeeze(stats(4,2,:)),squeeze(stats(4,3,:)), ...
    'markerLabel',td_longname, 'markerLegend', 'on', ...
    'styleSTD', '-', 'colOBS','k', 'markerObs','o', ...
    'markerSize',6, 'tickRMS',0:0.5:1, ...
    'tickRMSangle', 115, 'showlabelsRMS', 'on', ...
    'titleRMS','on', 'titleOBS','Obs');
t4=title('Surface dissolved Fe in PFZ','FontSize',14);
t4.Position(2) = t4.Position(2) + 0.3;

subplot(2,3,5)
[hp, ht, axl] = taylor_diagram(squeeze(stats(5,1,:)),squeeze(stats(5,2,:)),squeeze(stats(5,3,:)), ...
    'markerLabel',td_longname, 'markerLegend', 'on', ...
    'styleSTD', '-', 'colOBS','k', 'markerObs','o', ...
    'markerSize',6, 'tickRMS',0:0.5:1, ...
    'tickRMSangle', 115, 'showlabelsRMS', 'on', ...
    'titleRMS','on', 'titleOBS','Obs');
t5=title('Surface dissolved Fe in AZ','FontSize',14);
t5.Position(2) = t5.Position(2) + 0.3;

print(figs4,'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Figure\figs4','-djpeg','-r300');
