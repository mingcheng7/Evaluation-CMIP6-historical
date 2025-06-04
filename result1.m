% Evaluating the performance of CMIP6 models in simulating Southern Ocean
% biogeochemistry
% Author: Cheng, M., Ellwood, M., and Maher, N. 
% Last modified: 26/05/2025
% Chlorophyll related analysis (Section 3.1 and 3.2)
clear;

% Load data
chlfile_obs = 'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Regrided_data\Observation\plankton_copernicus_final.nc';
chlfile_acc = 'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Regrided_data\CMIP\chl_ACCESS-ESM1-5_final.nc';
chlfile_can = 'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Regrided_data\CMIP\chl_CanESM5_final.nc';
chlfile_ces = 'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Regrided_data\CMIP\chl_CESM2_final.nc';
chlfile_cmc = 'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Regrided_data\CMIP\chl_CMCC-ESM2_final.nc';
chlfile_cnr = 'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Regrided_data\CMIP\chl_CNRM-ESM2-1_final.nc';
chlfile_gfd = 'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Regrided_data\CMIP\chl_GFDL-ESM4_final.nc';
chlfile_ips = 'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Regrided_data\CMIP\chl_IPSL-CM6A-LR_final.nc';
chlfile_mir = 'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Regrided_data\CMIP\chl_MIROC-ES2L_final.nc';
chlfile_mhm = 'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Regrided_data\CMIP\chl_MPI-ESM-1-2-HAM_final.nc';
chlfile_mhr = 'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Regrided_data\CMIP\chl_MPI-ESM1-2-HR_final.nc';
chlfile_mlr = 'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Regrided_data\CMIP\chl_MPI-ESM1-2-LR_final.nc';
chlfile_nlm = 'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Regrided_data\CMIP\chl_NorESM2-LM_final.nc';
chlfile_nmm = 'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Regrided_data\CMIP\chl_NorESM2-MM_final.nc';
chlfile_uke = 'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Regrided_data\CMIP\chl_UKESM1-0-LL_final.nc';

model_abbs = {'acc','can','ces','cmc','cnr','gfd','ips','mir','mhm','mhr','mlr','nlm','nmm','uke'};
model_longname = {'ACCESS-ESM1-5','CanESM5','CESM2','CMCC-ESM2','CNRM-ESM2-1','GFDL-ESM4',...
    'IPSL-CM6A-LR','MIROC-ES2L','MPI-ESM-1-2-HAM','MPI-ESM1-2-HR',...
    'MPI-ESM1-2-LR','NorESM2-LM','NorESM2-MM','UKESM1-0-LL'};
lev_name = {'lev','lev','lev_partial','lev','lev','lev','olevel','lev','lev','lev','lev','lev','lev','lev'};

lon = ncread(chlfile_obs,'lon');
lat = ncread(chlfile_obs,'lat');
lev = ncread(chlfile_obs,'depth');
tim = ncread(chlfile_obs,'time');
chl_obs_read = ncread(chlfile_obs,'chl');
lonr = [lon;lon(1)];

[lat2,lon2] = meshgrid(lat,lon);
[lat2r,lon2r] = meshgrid(lat,lonr);

for i = 1:length(model_abbs)
    filename = eval(['chlfile_' model_abbs{i}]);
    chl_read = ncread(filename,'chl');
    lev_read = ncread(filename,lev_name{i});
    eval(['chl_' model_abbs{i} '_read=chl_read;']);
    eval(['lev_' model_abbs{i} '=lev_read;']);
end

% Optimise some depth mannually
lev = double(lev(2:end));
lev_mhm(1) = 5;
lev_mhr(1) = 5;
lev_mlr(1) = 5;

% Interplate CMIP6 data to Copernicus data grid size
chl_obs = chl_obs_read([1:360,1],:,2:end,:);
for i = 1:length(model_abbs)
    chl_regrid = NaN(size(chl_obs_read,1),size(chl_obs_read,2),size(chl_obs_read,3)-1,size(chl_obs_read,4));
    chl_read = eval(['chl_' model_abbs{i} '_read']);
    lev_read = eval(['lev_' model_abbs{i}]);
    for ji = 1:size(chl_regrid,1)
        for jj = 1:size(chl_regrid,2)
            for jk = 1:size(chl_regrid,4)
                chl_regrid(ji,jj,:,jk) = interp1(lev_read,squeeze(chl_read(ji,jj,:,jk)),lev)*1e6;
            end
        end
    end
    chl361 = chl_regrid([1:360,1],:,:,:);
    eval(['chl_' model_abbs{i} '=chl361;']);
end

chl_ips = chl_ips*1e-3;
chl_cnr = chl_cnr*1e-3;

% Fill blank for CMCC-ESM2, CanESM5, IPSL-CM6A-LR, MPI-ESM-1-2-HAM
chl_can_regrid = chl_can;
chl_cmc_regrid = chl_cmc;
chl_cnr_regrid = chl_cnr;
chl_ips_regrid = chl_ips;
chl_mhm_regrid = chl_mhm;
chl_mlr_regrid = chl_mlr;

chl_can(74,:,:,:) = mean(chl_can_regrid([73,75],:,:,:));
chl_cmc(73,:,:,:) = mean(chl_cmc_regrid([72,75],:,:,:));
chl_cmc(74,:,:,:) = mean(chl_cmc_regrid([72,75],:,:,:));
chl_cnr(73,:,:,:) = mean(chl_cnr_regrid([72,75],:,:,:));
chl_cnr(74,:,:,:) = mean(chl_cnr_regrid([72,75],:,:,:));
chl_ips(73,:,:,:) = mean(chl_ips_regrid([72,75],:,:,:));
chl_ips(74,:,:,:) = mean(chl_ips_regrid([72,75],:,:,:));

lon_fill_mhm = [147;148;148;149;149;150;150;151;151;152;152;153;153;153;...
    154;154;154;155;155;155;155;156;156;156;156;157;157;157;157;157];
lat_fill_mhm = [24;25;26;27;28;29;30;31;32;34;35;37;38;39;41;42;43;45;46;...
    47;48;50;51;52;53;56;57;58;59;60];

for f = 1:length(lon_fill_mhm)
    for ji = 1:length(lon_fill_mhm)
        for jj = 1:size(chl_mhm,3)
            for jk = 1:size(chl_mhm,4)
                chl_mhm(lon_fill_mhm(ji),lat_fill_mhm(ji),jj,jk) = ...
                    mean(chl_mhm_regrid([lon_fill_mhm(ji)-1,lon_fill_mhm(ji)+1],lat_fill_mhm(ji),jj,jk));
            end
        end
    end
end

for f = 1:length(lon_fill_mhm)
    for ji = 1:length(lon_fill_mhm)
        for jj = 1:size(chl_mlr,3)
            for jk = 1:size(chl_mlr,4)
                chl_mlr(lon_fill_mhm(ji),lat_fill_mhm(ji),jj,jk) = ...
                    mean(chl_mlr_regrid([lon_fill_mhm(ji)-1,lon_fill_mhm(ji)+1],lat_fill_mhm(ji),jj,jk));
            end
        end
    end
end

% Calculate mean chlorophyll in austral summer (DJJ) from 2000 to 2014
idx = ([1,2,12].*ones(15,1)+(0:12:14*12)')';
idx1 = idx(:);
chl_obs_djj_mean = nanmean(chl_obs(:,:,:,idx1),4);

for i = 1:length(model_abbs)
    chl_regrid = eval(['chl_' model_abbs{i}]);
    chl_djj_mean = mean(chl_regrid(:,:,:,idx1),4);
    chl_diff = chl_djj_mean - chl_obs_djj_mean;
    eval(['chl_' model_abbs{i} '_djj_mean=chl_djj_mean;']);
    eval(['chl_' model_abbs{i} '_diff=chl_diff;']);
end

% Express chl as a special scale for plotting
for i = 1:length(model_abbs)
    data = chl_obs_djj_mean;
    data(data>0.01 & data<=0.1) = (log10(data(data>0.01 & data<=0.1))+1)*0.2;
    data(data>0.1 & data<=0.2) = (data(data>0.1 & data<=0.2)-0.1)*2;
    data(data>1 & data<=1.5) = (data(data>1 & data<=1.5)-1)*0.4+1;
    data(data>1.5 & data<=2) = (data(data>1.5 & data<=2)-1.5)*0.4+1.2;
    data(data>2 & data<=5) = (data(data>2 & data<=5)-2)*0.0667+1.4;
    data(data>5 & data<=10) = (data(data>5 & data<=10)-1)*0.04+1.6;
    chl_obs_djj_mean_scale = data;
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

% Plot surface chlorophyll (Fig.1)
fig1=figure(1);
set(gcf,'Position',[100 100 1200 600])
sp0=subplot(3,6,1);
pos = get(gca, 'Position'); 
set(gca, 'Position', [pos(1) pos(2) pos(3)*1.15 pos(4)*1.15]); 
m_proj('stereographic','lat',-90,'long',0,'radius',60);
m_pcolor(lon2r,lat2r,chl_obs_djj_mean_scale(:,:,1))
shading flat;
colormap(sp0,chl_colmap);
m_line(line_stf(:,1),line_stf(:,2),'color','k','LineStyle','--','linewi',.5)
m_line(line_saf(:,1),line_saf(:,2),'color','k','LineStyle','--','linewi',.5)
m_line(line_pf(:,1),line_pf(:,2),'color','k','LineStyle','--','linewi',.5)
m_coast('patch','w');
m_coast('linewidth',1,'color','k');
m_grid('xtick',6,'xticklabel',[],'tickdir','out','ytick',[],'linest','none','xaxisloc','top','yaxislocation','middle','backcolor',[0.7 0.7 0.7]);
title('Copernicus','FontSize',12)
clim([-0.2 1.8])
c=colorbar;
c.Ticks = -0.2:0.2:1.8;
c.TickLabels = [0.01, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0, 1.5, 2, 5, 10]; 
c.Label.String = 'Chlorophyll (mg m^{-3})';
c.Label.FontSize = 11;
set(c,'Location','southoutside','Position',[0.14, 0.65, 0.1, 0.02])

for s = 1:length(splist)
    plotdata = eval(['chl_' model_abbs{s} '_diff']);
    sp=subplot(3,6,splist(s));
    pos = get(gca, 'Position'); 
    set(gca, 'Position', [pos(1) pos(2) pos(3)*1.15 pos(4)*1.15]); 
    m_proj('stereographic','lat',-90,'long',0,'radius',60);
    m_pcolor(lon2r,lat2r,plotdata(:,:,1))
    shading flat;
    colormap(sp,b2r(-2,2));
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
c1.Label.String = 'Chlorophyll (mg m^{-3})';
c1.Label.FontSize = 11;

print(fig1,'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Figure\fig1','-djpeg','-r600');

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

mbe_chl_os = nan(length(model_abbs),5);
for i = 1:length(model_abbs)
    chl_diff_os = eval(['chl_' model_abbs{i} '_diff(:,:,1)']);
    mbe_chl_os(i,1) = sum(chl_diff_os(~isnan(chl_diff_os)).*cosd(lat2r(~isnan(chl_diff_os))))...
        ./sum((cosd(lat2r(~isnan(chl_diff_os)))));
    mbe_chl_os(i,2) = sum(chl_diff_os(~isnan(chl_diff_os)&idx_stz).*cosd(lat2r(~isnan(chl_diff_os)&idx_stz)))...
        ./sum(cosd(lat2r(~isnan(chl_diff_os)&idx_stz)));
    mbe_chl_os(i,3) = sum(chl_diff_os(~isnan(chl_diff_os)&idx_saz).*cosd(lat2r(~isnan(chl_diff_os)&idx_saz)))...
        ./sum(cosd(lat2r(~isnan(chl_diff_os)&idx_saz)));
    mbe_chl_os(i,4) = sum(chl_diff_os(~isnan(chl_diff_os)&idx_pfz).*cosd(lat2r(~isnan(chl_diff_os)&idx_pfz)))...
        ./sum(cosd(lat2r(~isnan(chl_diff_os)&idx_pfz)));
    mbe_chl_os(i,5) = sum(chl_diff_os(~isnan(chl_diff_os)&idx_az).*cosd(lat2r(~isnan(chl_diff_os)&idx_az)))...
        ./sum(cosd(lat2r(~isnan(chl_diff_os)&idx_az)));
end

% Plot MBE (Fig2)
col = [0.4 0 0.6; 1 0.7 0.7; 0 0.8 0; 1 0.8 0; 0.4 0.7 1];
fig2=figure(2);
set(gcf,'Position',[100,100,1200,500])
b=bar(mbe_chl_os);
for c = 1:length(b)
    b(c).FaceColor = col(c,:);
end
ylabel('Mean Bias Error (mg m^{-3})')
set(gca,'XTickLabel',model_longname,'FontSize',12)
legend({'SO','STZ','SAZ','PFZ','AZ'})

print(fig2,'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Figure\fig2','-djpeg','-r600');

% Taylor Diagram
% Calculation of different parameters
tdlist = {'obs','acc','can','ces','cmc','cnr','gfd','ips','mir','mhm','mhr','mlr','nlm','nmm','uke'};
for i = 1:length(tdlist)
    chl_so = eval(['chl_' tdlist{i} '_djj_mean(:,:,1)']);
    chl_stz = chl_so(idx_stz);
    chl_saz = chl_so(idx_saz);
    chl_pfz = chl_so(idx_pfz);
    chl_az = chl_so(idx_az);
    eval(['chl_so_' tdlist{i} '=chl_so(:);']);
    eval(['chl_stz_' tdlist{i} '=chl_stz;']);
    eval(['chl_saz_' tdlist{i} '=chl_saz;']);
    eval(['chl_pfz_' tdlist{i} '=chl_pfz;']);
    eval(['chl_az_' tdlist{i} '=chl_az;']);
end

rglist = {'so','stz','saz','pfz','az'};

valid_idx_so = ~isnan(chl_so_obs) & ~isnan(chl_so_acc)& ~isnan(chl_so_can)  & ~isnan(chl_so_ces) & ...
    ~isnan(chl_so_cmc) & ~isnan(chl_so_cnr) & ~isnan(chl_so_gfd) & ~isnan(chl_so_ips) & ...
    ~isnan(chl_so_mir) & ~isnan(chl_so_mhm) & ~isnan(chl_so_mhr) & ~isnan(chl_so_mlr) & ...
    ~isnan(chl_so_nlm) & ~isnan(chl_so_nmm) & ~isnan(chl_so_uke);

valid_idx_stz = ~isnan(chl_stz_obs) & ~isnan(chl_stz_acc)& ~isnan(chl_stz_can)  & ~isnan(chl_stz_ces) & ...
    ~isnan(chl_stz_cmc) & ~isnan(chl_stz_cnr) & ~isnan(chl_stz_gfd) & ~isnan(chl_stz_ips) & ...
    ~isnan(chl_stz_mir) & ~isnan(chl_stz_mhm) & ~isnan(chl_stz_mhr) & ~isnan(chl_stz_mlr) & ...
    ~isnan(chl_stz_nlm) & ~isnan(chl_stz_nmm) & ~isnan(chl_stz_uke);

valid_idx_saz = ~isnan(chl_saz_obs) & ~isnan(chl_saz_acc)& ~isnan(chl_saz_can)  & ~isnan(chl_saz_ces) & ...
    ~isnan(chl_saz_cmc) & ~isnan(chl_saz_cnr) & ~isnan(chl_saz_gfd) & ~isnan(chl_saz_ips) & ...
    ~isnan(chl_saz_mir) & ~isnan(chl_saz_mhm) & ~isnan(chl_saz_mhr) & ~isnan(chl_saz_mlr) & ...
    ~isnan(chl_saz_nlm) & ~isnan(chl_saz_nmm) & ~isnan(chl_saz_uke);

valid_idx_pfz = ~isnan(chl_pfz_obs) & ~isnan(chl_pfz_acc)& ~isnan(chl_pfz_can)  & ~isnan(chl_pfz_ces) & ...
    ~isnan(chl_pfz_cmc) & ~isnan(chl_pfz_cnr) & ~isnan(chl_pfz_gfd) & ~isnan(chl_pfz_ips) & ...
    ~isnan(chl_pfz_mir) & ~isnan(chl_pfz_mhm) & ~isnan(chl_pfz_mhr) & ~isnan(chl_pfz_mlr) & ...
    ~isnan(chl_pfz_nlm) & ~isnan(chl_pfz_nmm) & ~isnan(chl_pfz_uke);

valid_idx_az = ~isnan(chl_az_obs) & ~isnan(chl_az_acc)& ~isnan(chl_az_can)  & ~isnan(chl_az_ces) & ...
    ~isnan(chl_az_cmc) & ~isnan(chl_az_cnr) & ~isnan(chl_az_gfd) & ~isnan(chl_az_ips) & ...
    ~isnan(chl_az_mir) & ~isnan(chl_az_mhm) & ~isnan(chl_az_mhr) & ~isnan(chl_az_mlr) & ...
    ~isnan(chl_az_nlm) & ~isnan(chl_az_nmm) & ~isnan(chl_az_uke);

for j = 1:length(rglist)
    for i = 1:length(tdlist)
        eval(['chl_' rglist{j} '_' tdlist{i} '_nonan=chl_' rglist{j} '_' tdlist{i} '(valid_idx_' rglist{j} ');']);
    end
end

stats = [];
for j = 1:length(rglist)
    for i = 1:length(model_abbs)
        S = SStats(eval(['chl_' rglist{j} '_obs_nonan']),eval(['chl_' rglist{j} '_' model_abbs{i} '_nonan']));
        stats(j,:,i+1) = S;
    end
    stats(j,:,1) = SStats(eval(['chl_' rglist{j} '_obs_nonan']),eval(['chl_' rglist{j} '_obs_nonan']));

    stdev = stats(j,2,1); 

    stats(j,2,:) = stats(j,2,:)/stdev; 
    stats(j,3,:) = stats(j,3,:)/stdev;
end

stats(:,1,:) = [];

td_longname = {' ','ACCESS-ESM1-5','CanESM5','CESM2','CMCC-ESM2','CNRM-ESM2-1','GFDL-ESM4',...
    'IPSL-CM6A-LR','MIROC-ES2L','MPI-ESM-1-2-HAM','MPI-ESM1-2-HR',...
    'MPI-ESM1-2-LR','NorESM2-LM','NorESM2-MM','UKESM1-0-LL'};

% Taylor Diagram for surface chlorophyll (Fig.S1)
figs1=figure(3);
set(gcf,'Position',[100 100 1500 1000])
subplot(2,3,1)
[hp, ht, axl] = taylor_diagram(squeeze(stats(1,1,:)),squeeze(stats(1,2,:)),squeeze(stats(1,3,:)), ...
    'markerLabel',td_longname,'markerLegend', 'on', ...
    'styleSTD', '-', 'colOBS','k', 'markerObs','o', ...
    'markerSize',6, 'tickRMS',0:2:6, ...
    'tickRMSangle', 115, 'showlabelsRMS', 'on', ...
    'titleRMS','on', 'titleOBS','Obs');
t1=title('Surface chlorophyll in SO','FontSize',14);
t1.Position(2) = t1.Position(2) + 0.8;

subplot(2,3,2)
[hp, ht, axl] = taylor_diagram(squeeze(stats(2,1,:)),squeeze(stats(2,2,:)),squeeze(stats(2,3,:)), ...
    'markerLabel',td_longname, 'markerLegend', 'on', ...
    'styleSTD', '-', 'colOBS','k', 'markerObs','o', ...
    'markerSize',6, 'tickRMS',0:1:2, ...
    'tickRMSangle', 115, 'showlabelsRMS', 'on', ...
    'titleRMS','on', 'titleOBS','Obs');
t2=title('Surface chlorophyll in STZ','FontSize',14);
t2.Position(2) = t2.Position(2) + 0.2;

subplot(2,3,3)
[hp, ht, axl] = taylor_diagram(squeeze(stats(3,1,:)),squeeze(stats(3,2,:)),squeeze(stats(3,3,:)), ...
    'markerLabel',td_longname, 'markerLegend', 'on', ...
    'styleSTD', '-', 'colOBS','k', 'markerObs','o', ...
    'markerSize',6, 'tickRMS',0:2:6, ...
    'tickRMSangle', 115, 'showlabelsRMS', 'on', ...
    'titleRMS','on', 'titleOBS','Obs');
t3=title('Surface chlorophyll in SAZ','FontSize',14);
t3.Position(2) = t3.Position(2) + 3;

subplot(2,3,4)
[hp, ht, axl] = taylor_diagram(squeeze(stats(4,1,:)),squeeze(stats(4,2,:)),squeeze(stats(4,3,:)), ...
    'markerLabel',td_longname, 'markerLegend', 'on', ...
    'styleSTD', '-', 'colOBS','k', 'markerObs','o', ...
    'markerSize',6, 'tickRMS',0:2:4, ...
    'tickRMSangle', 115, 'showlabelsRMS', 'on', ...
    'titleRMS','on', 'titleOBS','Obs');
t4=title('Surface chlorophyll in PFZ','FontSize',14);
t4.Position(2) = t4.Position(2) + 2;

subplot(2,3,5)
[hp, ht, axl] = taylor_diagram(squeeze(stats(5,1,:)),squeeze(stats(5,2,:)),squeeze(stats(5,3,:)), ...
    'markerLabel',td_longname, 'markerLegend', 'on', ...
    'styleSTD', '-', 'colOBS','k', 'markerObs','o', ...
    'markerSize',6, 'tickRMS',0:2:6, ...
    'tickRMSangle', 115, 'showlabelsRMS', 'on', ...
    'titleRMS','on', 'titleOBS','Obs');
t5=title('Surface chlorophyll in AZ','FontSize',14);
t5.Position(2) = t5.Position(2) + 3;

print(figs1,'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Figure\figs1','-djpeg','-r600');

% DCM plotting
% DCM analysis
dcmlist = {'obs','acc','can','ces','cmc','cnr','gfd','ips','mir','mhm','mhr','mlr','nlm','nmm','uke'};
dcmlist_longname = {'Copernicus','ACCESS-ESM1-5','CanESM5','CESM2','CMCC-ESM2','CNRM-ESM2-1','GFDL-ESM4',...
    'IPSL-CM6A-LR','MIROC-ES2L','MPI-ESM-1-2-HAM','MPI-ESM1-2-HR',...
    'MPI-ESM1-2-LR','NorESM2-LM','NorESM2-MM','UKESM1-0-LL'};
lev_obs = lev;

% Calculate mean chlorophyll in austral summer (DJJ) from 2000 to 2014
chl_obs_read_djj_mean = chl_obs_djj_mean;

for i = 1:length(model_abbs)
    chl_read = eval(['chl_' model_abbs{i} '_read']);
    chl_read_djj_mean = mean(chl_read([1:360,1],:,:,idx1),4)*1e6;
    eval(['chl_' model_abbs{i} '_read_djj_mean=chl_read_djj_mean;']);
end

chl_cnr_read_djj_mean = chl_cnr_read_djj_mean*1e-3;
chl_ips_read_djj_mean = chl_ips_read_djj_mean*1e-3;

% Fill blank for CMCC-ESM2, CanESM5, IPSL-CM6A-LR, MPI-ESM-1-2-HAM
chl_can_read_djj_mean_bac = chl_can_read_djj_mean;
chl_cmc_read_djj_mean_bac = chl_cmc_read_djj_mean;
chl_cnr_read_djj_mean_bac = chl_cnr_read_djj_mean;
chl_ips_read_djj_mean_bac = chl_ips_read_djj_mean;
chl_mhm_read_djj_mean_bac = chl_mhm_read_djj_mean;
chl_mlr_read_djj_mean_bac = chl_mlr_read_djj_mean;

chl_can_read_djj_mean(74,:,:,:) = mean(chl_can_read_djj_mean_bac([73,75],:,:,:));
chl_cmc_read_djj_mean(73,:,:,:) = mean(chl_cmc_read_djj_mean_bac([72,75],:,:,:));
chl_cmc_read_djj_mean(74,:,:,:) = mean(chl_cmc_read_djj_mean_bac([72,75],:,:,:));
chl_cnr_read_djj_mean(73,:,:,:) = mean(chl_cnr_read_djj_mean_bac([72,75],:,:,:));
chl_cnr_read_djj_mean(74,:,:,:) = mean(chl_cnr_read_djj_mean_bac([72,75],:,:,:));
chl_ips_read_djj_mean(73,:,:,:) = mean(chl_ips_read_djj_mean_bac([72,75],:,:,:));
chl_ips_read_djj_mean(74,:,:,:) = mean(chl_ips_read_djj_mean_bac([72,75],:,:,:));

lon_fill_mhm = [147;148;148;149;149;150;150;151;151;152;152;153;153;153;...
    154;154;154;155;155;155;155;156;156;156;156;157;157;157;157;157];
lat_fill_mhm = [24;25;26;27;28;29;30;31;32;34;35;37;38;39;41;42;43;45;46;...
    47;48;50;51;52;53;56;57;58;59;60];

for f = 1:length(lon_fill_mhm)
    for ji = 1:length(lon_fill_mhm)
        for jj = 1:size(chl_mhm_read_djj_mean,3)
            for jk = 1:size(chl_mhm_read_djj_mean,4)
                chl_mhm_read_djj_mean(lon_fill_mhm(ji),lat_fill_mhm(ji),jj,jk) = ...
                    mean(chl_mhm_read_djj_mean_bac([lon_fill_mhm(ji)-1,lon_fill_mhm(ji)+1],lat_fill_mhm(ji),jj,jk));
            end
        end
    end
end

for f = 1:length(lon_fill_mhm)
    for ji = 1:length(lon_fill_mhm)
        for jj = 1:size(chl_mlr_read_djj_mean,3)
            for jk = 1:size(chl_mlr_read_djj_mean,4)
                chl_mlr_read_djj_mean(lon_fill_mhm(ji),lat_fill_mhm(ji),jj,jk) = ...
                    mean(chl_mlr_read_djj_mean_bac([lon_fill_mhm(ji)-1,lon_fill_mhm(ji)+1],lat_fill_mhm(ji),jj,jk));
            end
        end
    end
end

for i = 1:length(dcmlist)
    data = eval(['chl_' dcmlist{i} '_read_djj_mean']);
    depth = eval(['lev_' dcmlist{i}]);
    [CHLM,CHLMind] = max(data,[],3); % Calculate the max
    CHLMind(isnan(CHLM)) = NaN;
    CHLM(CHLM<=data(:,:,1)*1.1) = 0;
    CHLMind(CHLM<=data(:,:,1)*1.1) = NaN;
    CHLM = squeeze(CHLM);
    CHLMind = squeeze(CHLMind);
    % chlorophyll concentration at each point and the depth index 
    DCM = CHLM; % Calculate the deep chlorophyll maxima, in other words, 
    DCMind = CHLMind;
    DCMd = NaN(size(DCM));

    % Calculate DCMs depth 
    for ji = 1:size(DCMind,1)
        for jj = 1:size(DCMind,2)
            if isnan(DCMind(ji,jj))
                DCMd(ji,jj) = NaN;
            elseif DCMind(ji,jj) == 1
                DCMd(ji,jj) = 0;
            else
                DCMd(ji,jj) = depth(DCMind(ji,jj));
            end
        end
    end

    % Export the variables
    eval(['dcm_' dcmlist{i} '=DCM;']);
    eval(['cmd_' dcmlist{i} '=DCMd;']);
end

% Express chl as a special scale for plotting
for i = 1:length(dcmlist)
    data = eval(['dcm_' dcmlist{i}]);
    data(data==0) = -0.2;
    data(data>0 & data<=0.01) = -0.1922;
    data(data>0.01 & data<=0.1) = (log10(data(data>0.01 & data<=0.1))+1)*0.2;
    data(data>0.1 & data<=0.2) = (data(data>0.1 & data<=0.2)-0.1)*2;
    data(data>1 & data<=1.5) = (data(data>1 & data<=1.5)-1)*0.4+1;
    data(data>1.5 & data<=2) = (data(data>1.5 & data<=2)-1.5)*0.4+1.2;
    data(data>2 & data<=5) = (data(data>2 & data<=5)-2)*0.0667+1.4;
    data(data>5 & data<=10) = (data(data>5 & data<=10)-1)*0.04+1.6;
    eval(['dcm_' dcmlist{i} '_scale=data;']);
end

splist2 = [1;2;3;4;5;7;8;9;10;11;13;14;15;16;17];
chl_nan = chl_colmap;
chl_nan(1,:) = 1;

% Plot DCM (Fig.9)
fig9=figure(4);
set(gcf,'Position',[100 100 1200 600])
for s = 1:length(splist2)
    plotdata = eval(['dcm_' dcmlist{s} '_scale']);
    sp=subplot(3,6,splist2(s));
    pos = get(gca, 'Position'); 
    set(gca, 'Position', [pos(1) pos(2) pos(3)*1.15 pos(4)*1.15]); 
    m_proj('stereographic','lat',-90,'long',0,'radius',60);
    m_pcolor(lon2r,lat2r,plotdata)
    shading flat;
    colormap(sp,chl_nan);
    m_line(line_stf(:,1),line_stf(:,2),'color','k','LineStyle','--','linewi',.5)
    m_line(line_saf(:,1),line_saf(:,2),'color','k','LineStyle','--','linewi',.5)
    m_line(line_pf(:,1),line_pf(:,2),'color','k','LineStyle','--','linewi',.5)
    m_coast('patch','w');
    m_coast('linewidth',1,'color','k');
    m_grid('xtick',6,'xticklabel',[],'tickdir','out','ytick',[],'linest','none','xaxisloc','top','yaxislocation','middle','backcolor',[0.7 0.7 0.7]);
    title(dcmlist_longname{s},'FontSize',12)
    clim([-0.2 1.8])
end

c=colorbar;
c.Ticks = -0.2:0.2:1.8;
c.TickLabels = [0.01, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0, 1.5, 2, 5, 10]; 
c.Label.String = 'Chlorophyll (mg m^{-3})';
c.Label.FontSize = 11;
set(c,'Location','eastoutside','Position',[0.8, 0.4, 0.01, 0.3])

print(fig9,'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Figure\fig9','-djpeg','-r600');

% Plot the number of grids
dcmgrid = nan(length(dcmlist),5);
for i = 1:length(dcmlist)
    dcm = eval(['dcm_' dcmlist{i}]);
    dcmgrid(i,1) = sum((dcm>0).*cosd(lat2r),'all')...
        ./sum((cosd(lat2r(~isnan(dcm)))));
    dcmgrid(i,2) = sum((dcm>0&idx_stz).*cosd(lat2r),'all')...
        ./sum(cosd(lat2r(~isnan(dcm)&idx_stz)));
    dcmgrid(i,3) = sum((dcm>0&idx_saz).*cosd(lat2r),'all')...
        ./sum(cosd(lat2r(~isnan(dcm)&idx_saz)));
    dcmgrid(i,4) = sum((dcm>0&idx_pfz).*cosd(lat2r),'all')...
        ./sum(cosd(lat2r(~isnan(dcm)&idx_pfz)));
    dcmgrid(i,5) = sum((dcm>0&idx_az).*cosd(lat2r),'all')...
        ./sum(cosd(lat2r(~isnan(dcm)&idx_az)));
end

% Plot (Fig.10)
col = [0.4 0 0.6; 1 0.7 0.7; 0 0.8 0; 1 0.8 0; 0.4 0.7 1];
fig10=figure(5);
set(gcf,'Position',[100,100,1200,500])
b=bar(dcmgrid*100);
for c = 1:length(b)
    b(c).FaceColor = col(c,:);
end
ylim([0 100])
ylabel('Percentage (%)')
set(gca,'XTickLabel',dcmlist_longname,'FontSize',12)
legend({'SO','STZ','SAZ','PFZ','AZ'})

print(fig10,'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Figure\fig10','-djpeg','-r600');
