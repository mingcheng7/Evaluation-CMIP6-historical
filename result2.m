% Evaluating the performance of CMIP6 models in simulating Southern Ocean
% biogeochemistry
% Author: Cheng, M., Ellwood, M., and Maher, N. 
% Last modified: 26/05/2025
% Nitrate related analysis (Section 3.1)
clear;

% Load data
no3file_obs = 'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Regrided_data\Observation\no3_woa_final.nc';
no3file_acc = 'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Regrided_data\CMIP\no3_ACCESS-ESM1-5_final.nc';
no3file_can = 'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Regrided_data\CMIP\no3_CanESM5_final.nc';
no3file_ces = 'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Regrided_data\CMIP\no3_CESM2_final.nc';
no3file_cmc = 'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Regrided_data\CMIP\no3_CMCC-ESM2_final.nc';
no3file_cnr = 'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Regrided_data\CMIP\no3_CNRM-ESM2-1_final.nc';
no3file_gfd = 'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Regrided_data\CMIP\no3_GFDL-ESM4_final.nc';
no3file_ips = 'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Regrided_data\CMIP\no3_IPSL-CM6A-LR_final.nc';
no3file_mir = 'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Regrided_data\CMIP\no3_MIROC-ES2L_final.nc';
no3file_mhm = 'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Regrided_data\CMIP\no3_MPI-ESM-1-2-HAM_final.nc';
no3file_mhr = 'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Regrided_data\CMIP\no3_MPI-ESM1-2-HR_final.nc';
no3file_mlr = 'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Regrided_data\CMIP\no3_MPI-ESM1-2-LR_final.nc';
no3file_nlm = 'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Regrided_data\CMIP\no3_NorESM2-LM_final.nc';
no3file_nmm = 'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Regrided_data\CMIP\no3_NorESM2-MM_final.nc';
no3file_uke = 'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Regrided_data\CMIP\no3_UKESM1-0-LL_final.nc';

model_abbs = {'acc','can','ces','cmc','cnr','gfd','ips','mir','mhm','mhr','mlr','nlm','nmm','uke'};
model_longname = {'ACCESS-ESM1-5','CanESM5','CESM2','CMCC-ESM2','CNRM-ESM2-1','GFDL-ESM4',...
    'IPSL-CM6A-LR','MIROC-ES2L','MPI-ESM-1-2-HAM','MPI-ESM1-2-HR',...
    'MPI-ESM1-2-LR','NorESM2-LM','NorESM2-MM','UKESM1-0-LL'};
lev_name = {'lev','lev','lev','lev','lev','lev','olevel','lev','lev','lev','lev','lev','lev','lev'};

lon = ncread(no3file_obs,'lon');
lat = ncread(no3file_obs,'lat');
lev = ncread(no3file_obs,'depth');
tim = ncread(no3file_obs,'time');
no3_obs_read = ncread(no3file_obs,'n_an');

lonr = [lon;lon(1)];

[lat2,lon2] = meshgrid(lat,lon);
[lat2r,lon2r] = meshgrid(lat,lonr);

for i = 1:length(model_abbs)
    filename = eval(['no3file_' model_abbs{i}]);
    no3_read = ncread(filename,'no3');
    lev_read = ncread(filename,lev_name{i});
    eval(['no3_' model_abbs{i} '_read=no3_read;']);
    eval(['lev_' model_abbs{i} '=lev_read;']);
end

% Optimise some depth mannually
lev = double(lev(2:end));
lev_mhm(1) = 5;
lev_mhr(1) = 5;
lev_mlr(1) = 5;

% Interplate CMIP6 data to Copernicus data grid size
no3_obs = no3_obs_read([1:360,1],:,2:end,:);
for i = 1:length(model_abbs)
    no3_regrid = NaN(size(no3_obs,1)-1,size(no3_obs,2),size(no3_obs,3),180);
    no3_read = eval(['no3_' model_abbs{i} '_read']);
    lev_read = eval(['lev_' model_abbs{i}]);
    for ji = 1:size(no3_regrid,1)
        for jj = 1:size(no3_regrid,2)
            for jk = 1:size(no3_regrid,4)
                no3_regrid(ji,jj,:,jk) = interp1(lev_read,squeeze(no3_read(ji,jj,:,jk)),lev)*1000;
            end
        end
    end
    no3361 = no3_regrid([1:360,1],:,:,:);
    eval(['no3_' model_abbs{i} '=no3361;']);
end

% Fill blank for CMCC-ESM2, CanESM5, IPSL-CM6A-LR, MPI-ESM-1-2-HAM
no3_can_regrid = no3_can;
no3_cmc_regrid = no3_cmc;
no3_cnr_regrid = no3_cnr;
no3_ips_regrid = no3_ips;
no3_mhm_regrid = no3_mhm;
no3_mlr_regrid = no3_mlr;

no3_can(74,:,:,:) = mean(no3_can_regrid([73,75],:,:,:));
no3_cmc(73,:,:,:) = mean(no3_cmc_regrid([72,75],:,:,:));
no3_cmc(74,:,:,:) = mean(no3_cmc_regrid([72,75],:,:,:));
no3_cnr(73,:,:,:) = mean(no3_cnr_regrid([72,75],:,:,:));
no3_cnr(74,:,:,:) = mean(no3_cnr_regrid([72,75],:,:,:));
no3_ips(73,:,:,:) = mean(no3_ips_regrid([72,75],:,:,:));
no3_ips(74,:,:,:) = mean(no3_ips_regrid([72,75],:,:,:));

lon_fill_mhm = [147;148;148;149;149;150;150;151;151;152;152;153;153;153;...
    154;154;154;155;155;155;155;156;156;156;156;157;157;157;157;157];
lat_fill_mhm = [24;25;26;27;28;29;30;31;32;34;35;37;38;39;41;42;43;45;46;...
    47;48;50;51;52;53;56;57;58;59;60];

for f = 1:length(lon_fill_mhm)
    for ji = 1:length(lon_fill_mhm)
        for jj = 1:size(no3_mhm,3)
            for jk = 1:size(no3_mhm,4)
                no3_mhm(lon_fill_mhm(ji),lat_fill_mhm(ji),jj,jk) = ...
                    mean(no3_mhm_regrid([lon_fill_mhm(ji)-1,lon_fill_mhm(ji)+1],lat_fill_mhm(ji),jj,jk));
            end
        end
    end
end

for f = 1:length(lon_fill_mhm)
    for ji = 1:length(lon_fill_mhm)
        for jj = 1:size(no3_mlr,3)
            for jk = 1:size(no3_mlr,4)
                no3_mlr(lon_fill_mhm(ji),lat_fill_mhm(ji),jj,jk) = ...
                    mean(no3_mlr_regrid([lon_fill_mhm(ji)-1,lon_fill_mhm(ji)+1],lat_fill_mhm(ji),jj,jk));
            end
        end
    end
end

% Calculate mean nitrate in austral summer (DJJ) from 2000 to 2014
idx = ([1,2,12].*ones(15,1)+(0:12:14*12)')';
idx1 = idx(:);
no3_obs_djj_mean = mean(no3_obs(:,:,:,[1;2;12]),4);

for i = 1:length(model_abbs)
    no3_regrid = eval(['no3_' model_abbs{i}]);
    no3_djj_mean = mean(no3_regrid(:,:,:,idx1),4);
    no3_diff = no3_djj_mean - no3_obs_djj_mean;
    eval(['no3_' model_abbs{i} '_djj_mean=no3_djj_mean;']);
    eval(['no3_' model_abbs{i} '_diff=no3_diff;']);
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

% Plot surface nitrate (Fig.3)
fig3=figure(1);
set(gcf,'Position',[100 100 1200 600])
sp0=subplot(3,6,1);
pos = get(gca, 'Position'); 
set(gca, 'Position', [pos(1) pos(2) pos(3)*1.15 pos(4)*1.15]); 
m_proj('stereographic','lat',-90,'long',0,'radius',60);
m_pcolor(lon2r,lat2r,no3_obs_djj_mean(:,:,1))
shading flat;
colormap(sp0,flipud(m_colmap('blues')));
m_line(line_stf(:,1),line_stf(:,2),'color','k','LineStyle','--','linewi',.5)
m_line(line_saf(:,1),line_saf(:,2),'color','k','LineStyle','--','linewi',.5)
m_line(line_pf(:,1),line_pf(:,2),'color','k','LineStyle','--','linewi',.5)
m_coast('patch','w');
m_coast('linewidth',1,'color','k');
m_grid('xtick',6,'xticklabel',[],'tickdir','out','ytick',[],'linest','none','xaxisloc','top','yaxislocation','middle','backcolor',[0.7 0.7 0.7]);
title('WOA','FontSize',12)
clim([0 30])
c=colorbar;
c.Label.String = 'Nitrate (mmol/m^{3})';
c.Label.FontSize = 11;
set(c,'Location','southoutside','Position',[0.14,0.65,0.1,0.02])

for s = 1:length(splist)
    plotdata = eval(['no3_' model_abbs{s} '_diff']);
    sp=subplot(3,6,splist(s));
    pos = get(gca, 'Position'); 
    set(gca, 'Position', [pos(1) pos(2) pos(3)*1.15 pos(4)*1.15]); 
    m_proj('stereographic','lat',-90,'long',0,'radius',60);
    m_pcolor(lon2r,lat2r,plotdata(:,:,1))
    shading flat;
    colormap(sp,b2r(-30,30));
    m_line(line_stf(:,1),line_stf(:,2),'color','k','LineStyle','--','linewi',.5)
    m_line(line_saf(:,1),line_saf(:,2),'color','k','LineStyle','--','linewi',.5)
    m_line(line_pf(:,1),line_pf(:,2),'color','k','LineStyle','--','linewi',.5)
    m_coast('patch','w');
    m_coast('linewidth',1,'color','k');
    m_grid('xtick',6,'xticklabel',[],'tickdir','out','ytick',[],'linest','none','xaxisloc','top','yaxislocation','middle','backcolor',[0.7 0.7 0.7]);
    title(model_longname{s},'FontSize',12)
end
c1=colorbar;
set(c1,'Location','eastoutside','Position',[0.82, 0.12, 0.01, 0.25])
c1.Label.String = 'Nitrate (mmol/m^{3})';
c1.Label.FontSize = 11;

print(fig3,'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Figure\fig3','-djpeg','-r300');

% Calculate mean bias error for 4 zones and total SO
lon2p = lon2;
lon2p(lon2p>180) = lon2p(lon2p>180) - 360;
so_boundary1 = [-180,-30;180,-30];
so_boundary2 = [-180,-90;180,-90];
poly_stz = [so_boundary1;flipud(line_stf)];
poly_saz = [line_stf;flipud(line_saf)];
poly_pfz = [line_saf;flipud(line_pf)];
poly_az = [line_pf;flipud(so_boundary2)];

idx_stz = inpolygon(lon2r,lat2r,poly_stz(:,1),poly_stz(:,2));
idx_saz = inpolygon(lon2r,lat2r,poly_saz(:,1),poly_saz(:,2));
idx_pfz = inpolygon(lon2r,lat2r,poly_pfz(:,1),poly_pfz(:,2));
idx_az = inpolygon(lon2r,lat2r,poly_az(:,1),poly_az(:,2));

mbe_no3_os = nan(length(model_abbs),5);
for i = 1:length(model_abbs)
    no3_diff_os = eval(['no3_' model_abbs{i} '_diff(:,:,1)']);
    mbe_no3_os(i,1) = sum(no3_diff_os(~isnan(no3_diff_os)).*cosd(lat2r(~isnan(no3_diff_os))))...
        ./sum((cosd(lat2r(~isnan(no3_diff_os)))));
    mbe_no3_os(i,2) = sum(no3_diff_os(~isnan(no3_diff_os)&idx_stz).*cosd(lat2r(~isnan(no3_diff_os)&idx_stz)))...
        ./sum(cosd(lat2r(~isnan(no3_diff_os)&idx_stz)));
    mbe_no3_os(i,3) = sum(no3_diff_os(~isnan(no3_diff_os)&idx_saz).*cosd(lat2r(~isnan(no3_diff_os)&idx_saz)))...
        ./sum(cosd(lat2r(~isnan(no3_diff_os)&idx_saz)));
    mbe_no3_os(i,4) = sum(no3_diff_os(~isnan(no3_diff_os)&idx_pfz).*cosd(lat2r(~isnan(no3_diff_os)&idx_pfz)))...
        ./sum(cosd(lat2r(~isnan(no3_diff_os)&idx_pfz)));
    mbe_no3_os(i,5) = sum(no3_diff_os(~isnan(no3_diff_os)&idx_az).*cosd(lat2r(~isnan(no3_diff_os)&idx_az)))...
        ./sum(cosd(lat2r(~isnan(no3_diff_os)&idx_az)));
end

% Plot MBE (Fig.4)
col = [0.4 0 0.6; 1 0.7 0.7; 0 0.8 0; 1 0.8 0; 0.4 0.7 1];
fig4=figure(2);
set(gcf,'Position',[100,100,1200,500])
b=bar(mbe_no3_os);
for c = 1:length(b)
    b(c).FaceColor = col(c,:);
end
ylabel('Mean Bias Error (mmol/m^{3})')
set(gca,'XTickLabel',model_longname,'FontSize',12)
legend({'SO','STZ','SAZ','PFZ','AZ'},'Location','southeast')

print(fig4,'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Figure\fig4','-djpeg','-r300');

% Taylor Diagram
% Calculation of different parameters
tdlist = {'obs','acc','can','ces','cmc','cnr','gfd','ips','mir','mhm','mhr','mlr','nlm','nmm','uke'};
for i = 1:length(tdlist)
    no3_so = eval(['no3_' tdlist{i} '_djj_mean(:,:,1)']);
    no3_stz = no3_so(idx_stz);
    no3_saz = no3_so(idx_saz);
    no3_pfz = no3_so(idx_pfz);
    no3_az = no3_so(idx_az);
    eval(['no3_so_' tdlist{i} '=no3_so(:);']);
    eval(['no3_stz_' tdlist{i} '=no3_stz;']);
    eval(['no3_saz_' tdlist{i} '=no3_saz;']);
    eval(['no3_pfz_' tdlist{i} '=no3_pfz;']);
    eval(['no3_az_' tdlist{i} '=no3_az;']);
end

rglist = {'so','stz','saz','pfz','az'};

valid_idx_so = ~isnan(no3_so_obs) & ~isnan(no3_so_acc)& ~isnan(no3_so_can)  & ~isnan(no3_so_ces) & ...
    ~isnan(no3_so_cmc) & ~isnan(no3_so_cnr) & ~isnan(no3_so_gfd) & ~isnan(no3_so_ips) & ...
    ~isnan(no3_so_mir) & ~isnan(no3_so_mhm) & ~isnan(no3_so_mhr) & ~isnan(no3_so_mlr) & ...
    ~isnan(no3_so_nlm) & ~isnan(no3_so_nmm) & ~isnan(no3_so_uke);

valid_idx_stz = ~isnan(no3_stz_obs) & ~isnan(no3_stz_acc)& ~isnan(no3_stz_can)  & ~isnan(no3_stz_ces) & ...
    ~isnan(no3_stz_cmc) & ~isnan(no3_stz_cnr) & ~isnan(no3_stz_gfd) & ~isnan(no3_stz_ips) & ...
    ~isnan(no3_stz_mir) & ~isnan(no3_stz_mhm) & ~isnan(no3_stz_mhr) & ~isnan(no3_stz_mlr) & ...
    ~isnan(no3_stz_nlm) & ~isnan(no3_stz_nmm) & ~isnan(no3_stz_uke);

valid_idx_saz = ~isnan(no3_saz_obs) & ~isnan(no3_saz_acc)& ~isnan(no3_saz_can)  & ~isnan(no3_saz_ces) & ...
    ~isnan(no3_saz_cmc) & ~isnan(no3_saz_cnr) & ~isnan(no3_saz_gfd) & ~isnan(no3_saz_ips) & ...
    ~isnan(no3_saz_mir) & ~isnan(no3_saz_mhm) & ~isnan(no3_saz_mhr) & ~isnan(no3_saz_mlr) & ...
    ~isnan(no3_saz_nlm) & ~isnan(no3_saz_nmm) & ~isnan(no3_saz_uke);

valid_idx_pfz = ~isnan(no3_pfz_obs) & ~isnan(no3_pfz_acc)& ~isnan(no3_pfz_can)  & ~isnan(no3_pfz_ces) & ...
    ~isnan(no3_pfz_cmc) & ~isnan(no3_pfz_cnr) & ~isnan(no3_pfz_gfd) & ~isnan(no3_pfz_ips) & ...
    ~isnan(no3_pfz_mir) & ~isnan(no3_pfz_mhm) & ~isnan(no3_pfz_mhr) & ~isnan(no3_pfz_mlr) & ...
    ~isnan(no3_pfz_nlm) & ~isnan(no3_pfz_nmm) & ~isnan(no3_pfz_uke);

valid_idx_az = ~isnan(no3_az_obs) & ~isnan(no3_az_acc)& ~isnan(no3_az_can)  & ~isnan(no3_az_ces) & ...
    ~isnan(no3_az_cmc) & ~isnan(no3_az_cnr) & ~isnan(no3_az_gfd) & ~isnan(no3_az_ips) & ...
    ~isnan(no3_az_mir) & ~isnan(no3_az_mhm) & ~isnan(no3_az_mhr) & ~isnan(no3_az_mlr) & ...
    ~isnan(no3_az_nlm) & ~isnan(no3_az_nmm) & ~isnan(no3_az_uke);

for j = 1:length(rglist)
    for i = 1:length(tdlist)
        eval(['no3_' rglist{j} '_' tdlist{i} '_nonan=no3_' rglist{j} '_' tdlist{i} '(valid_idx_' rglist{j} ');']);
    end
end

stats = [];
for j = 1:length(rglist)
    for i = 1:length(model_abbs)
        S = SStats(eval(['no3_' rglist{j} '_obs_nonan']),eval(['no3_' rglist{j} '_' model_abbs{i} '_nonan']));
        stats(j,:,i+1) = S;
    end
    stats(j,:,1) = SStats(eval(['no3_' rglist{j} '_obs_nonan']),eval(['no3_' rglist{j} '_obs_nonan']));

    stdev = stats(j,2,1); 

    stats(j,2,:) = stats(j,2,:)/stdev; 
    stats(j,3,:) = stats(j,3,:)/stdev;
end

stats(:,1,:) = [];

% Taylor Diagram for surface nitrate (Fig.S2)
td_longname = {' ','ACCESS-ESM1-5','CanESM5','CESM2','CMCC-ESM2','CNRM-ESM2-1','GFDL-ESM4',...
    'IPSL-CM6A-LR','MIROC-ES2L','MPI-ESM-1-2-HAM','MPI-ESM1-2-HR',...
    'MPI-ESM1-2-LR','NorESM2-LM','NorESM2-MM','UKESM1-0-LL'};

figs2=figure(3);
set(gcf,'Position',[100 100 1500 1000])
subplot(2,3,1)
[hp, ht, axl] = taylor_diagram(squeeze(stats(1,1,:)),squeeze(stats(1,2,:)),squeeze(stats(1,3,:)), ...
    'markerLabel',td_longname,'markerLegend', 'on', ...
    'styleSTD', '-', 'colOBS','k', 'markerObs','o', ...
    'markerSize',6, 'tickRMS',0:0.5:1, ...
    'tickRMSangle', 115, 'showlabelsRMS', 'on', ...
    'titleRMS','on', 'titleOBS','Obs');
t1=title('Surface nitrate in SO','FontSize',14);
t1.Position(2) = t1.Position(2) + 0.12;

subplot(2,3,2)
[hp, ht, axl] = taylor_diagram(squeeze(stats(2,1,:)),squeeze(stats(2,2,:)),squeeze(stats(2,3,:)), ...
    'markerLabel',td_longname, 'markerLegend', 'on', ...
    'styleSTD', '-', 'colOBS','k', 'markerObs','o', ...
    'markerSize',6, 'tickRMS',0:1:2, ...
    'tickRMSangle', 115, 'showlabelsRMS', 'on', ...
    'titleRMS','on', 'titleOBS','Obs');
t2=title('Surface nitrate in STZ','FontSize',14);
t2.Position(2) = t2.Position(2) + 0.24;

subplot(2,3,3)
[hp, ht, axl] = taylor_diagram(squeeze(stats(3,1,:)),squeeze(stats(3,2,:)),squeeze(stats(3,3,:)), ...
    'markerLabel',td_longname, 'markerLegend', 'on', ...
    'styleSTD', '-', 'colOBS','k', 'markerObs','o', ...
    'markerSize',6, 'tickRMS',0:0.5:1, ...
    'tickRMSangle', 115, 'showlabelsRMS', 'on', ...
    'titleRMS','on', 'titleOBS','Obs');
t3=title('Surface nitrate in SAZ','FontSize',14);
t3.Position(2) = t3.Position(2) + 0.08;

subplot(2,3,4)
[hp, ht, axl] = taylor_diagram(squeeze(stats(4,1,:)),squeeze(stats(4,2,:)),squeeze(stats(4,3,:)), ...
    'markerLabel',td_longname, 'markerLegend', 'on', ...
    'styleSTD', '-', 'colOBS','k', 'markerObs','o', ...
    'markerSize',6, 'tickRMS',0:0.5:1, ...
    'tickRMSangle', 115, 'showlabelsRMS', 'on', ...
    'titleRMS','on', 'titleOBS','Obs');
t4=title('Surface nitrate in PFZ','FontSize',14);
t4.Position(2) = t4.Position(2) + 0.12;

subplot(2,3,5)
[hp, ht, axl] = taylor_diagram(squeeze(stats(5,1,:)),squeeze(stats(5,2,:)),squeeze(stats(5,3,:)), ...
    'markerLabel',td_longname, 'markerLegend', 'on', ...
    'styleSTD', '-', 'colOBS','k', 'markerObs','o', ...
    'markerSize',6, 'tickRMS',0:1:3, ...
    'tickRMSangle', 115, 'showlabelsRMS', 'on', ...
    'titleRMS','on', 'titleOBS','Obs');
t5=title('Surface nitrate in AZ','FontSize',14);
t5.Position(2) = t5.Position(2) + 1;

print(figs2,'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Figure\figs2','-djpeg','-r300');
