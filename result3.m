% Evaluating the performance of CMIP6 models in simulating Southern Ocean
% biogeochemistry
% Author: Cheng, M., Ellwood, M., and Maher, N. 
% Last modified: 26/05/2025
% Silicate related analysis (Section 3.1)
clear;

% Load data
sifile_obs = 'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Regrided_data\Observation\si_woa_final.nc';
sifile_ces = 'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Regrided_data\CMIP\sil_CESM2_final.nc';
sifile_cmc = 'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Regrided_data\CMIP\sil_CMCC-ESM2_final.nc';
sifile_cnr = 'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Regrided_data\CMIP\sil_CNRM-ESM2-1_final.nc';
sifile_gfd = 'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Regrided_data\CMIP\sil_GFDL-ESM4_final.nc';
sifile_ips = 'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Regrided_data\CMIP\sil_IPSL-CM6A-LR_final.nc';
sifile_mhm = 'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Regrided_data\CMIP\sil_MPI-ESM-1-2-HAM_final.nc';
sifile_mhr = 'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Regrided_data\CMIP\sil_MPI-ESM1-2-HR_final.nc';
sifile_mlr = 'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Regrided_data\CMIP\sil_MPI-ESM1-2-LR_final.nc';
sifile_nlm = 'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Regrided_data\CMIP\sil_NorESM2-LM_final.nc';
sifile_nmm = 'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Regrided_data\CMIP\sil_NorESM2-MM_final.nc';
sifile_uke = 'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Regrided_data\CMIP\sil_UKESM1-0-LL_final.nc';

model_abbs = {'ces','cmc','cnr','gfd','ips','mhm','mhr','mlr','nlm','nmm','uke'};
model_longname = {'CESM2','CMCC-ESM2','CNRM-ESM2-1','GFDL-ESM4','IPSL-CM6A-LR','MPI-ESM-1-2-HAM',...
    'MPI-ESM1-2-HR','MPI-ESM1-2-LR','NorESM2-LM','NorESM2-MM','UKESM1-0-LL'};
lev_name = {'lev','lev','lev','lev','olevel','lev','lev','lev','lev','lev','lev'};

lon = ncread(sifile_obs,'lon');
lat = ncread(sifile_obs,'lat');
lev = ncread(sifile_obs,'depth');
tim = ncread(sifile_obs,'time');
si_obs_read = ncread(sifile_obs,'i_an');

lonr = [lon;lon(1)];

[lat2,lon2] = meshgrid(lat,lon);
[lat2r,lon2r] = meshgrid(lat,lonr);

for i = 1:length(model_abbs)
    filename = eval(['sifile_' model_abbs{i}]);
    si_read = ncread(filename,'si');
    lev_read = ncread(filename,lev_name{i});
    eval(['si_' model_abbs{i} '_read=si_read;']);
    eval(['lev_' model_abbs{i} '=lev_read;']);
end

% Optimise some depth mannually
lev = double(lev(2:end));
lev_mhm(1) = 5;
lev_mhr(1) = 5;
lev_mlr(1) = 5;

% Interplate CMIP6 data to Copernicus data grid size
si_obs = si_obs_read([1:360,1],:,2:end,:);
for i = 1:length(model_abbs)
    si_regrid = NaN(size(si_obs,1)-1,size(si_obs,2),size(si_obs,3),180);
    si_read = eval(['si_' model_abbs{i} '_read']);
    lev_read = eval(['lev_' model_abbs{i}]);
    for ji = 1:size(si_regrid,1)
        for jj = 1:size(si_regrid,2)
            for jk = 1:size(si_regrid,4)
                si_regrid(ji,jj,:,jk) = interp1(lev_read,squeeze(si_read(ji,jj,:,jk)),lev)*1000;
            end
        end
    end
    si361 = si_regrid([1:360,1],:,:,:);
    eval(['si_' model_abbs{i} '=si361;']);
end

% Fill blank for CMCC-ESM2, CanESM5, IPSL-CM6A-LR, MPI-ESM-1-2-HAM
si_cmc_regrid = si_cmc;
si_cnr_regrid = si_cnr;
si_ips_regrid = si_ips;
si_mhm_regrid = si_mhm;
si_mlr_regrid = si_mlr;

si_cmc(73,:,:,:) = mean(si_cmc_regrid([72,75],:,:,:));
si_cmc(74,:,:,:) = mean(si_cmc_regrid([72,75],:,:,:));
si_cnr(73,:,:,:) = mean(si_cnr_regrid([72,75],:,:,:));
si_cnr(74,:,:,:) = mean(si_cnr_regrid([72,75],:,:,:));
si_ips(73,:,:,:) = mean(si_ips_regrid([72,75],:,:,:));
si_ips(74,:,:,:) = mean(si_ips_regrid([72,75],:,:,:));

lon_fill_mhm = [147;148;148;149;149;150;150;151;151;152;152;153;153;153;...
    154;154;154;155;155;155;155;156;156;156;156;157;157;157;157;157];
lat_fill_mhm = [24;25;26;27;28;29;30;31;32;34;35;37;38;39;41;42;43;45;46;...
    47;48;50;51;52;53;56;57;58;59;60];

for f = 1:length(lon_fill_mhm)
    for ji = 1:length(lon_fill_mhm)
        for jj = 1:size(si_mhm,3)
            for jk = 1:size(si_mhm,4)
                si_mhm(lon_fill_mhm(ji),lat_fill_mhm(ji),jj,jk) = ...
                    mean(si_mhm_regrid([lon_fill_mhm(ji)-1,lon_fill_mhm(ji)+1],lat_fill_mhm(ji),jj,jk));
            end
        end
    end
end

for f = 1:length(lon_fill_mhm)
    for ji = 1:length(lon_fill_mhm)
        for jj = 1:size(si_mlr,3)
            for jk = 1:size(si_mlr,4)
                si_mlr(lon_fill_mhm(ji),lat_fill_mhm(ji),jj,jk) = ...
                    mean(si_mlr_regrid([lon_fill_mhm(ji)-1,lon_fill_mhm(ji)+1],lat_fill_mhm(ji),jj,jk));
            end
        end
    end
end

% Calculate mean Si in austral summer (DJJ) from 2000 to 2014
idx = ([1,2,12].*ones(15,1)+(0:12:14*12)')';
idx1 = idx(:);
si_obs_djj_mean = mean(si_obs(:,:,:,[1;2;12]),4);

for i = 1:length(model_abbs)
    si_regrid = eval(['si_' model_abbs{i}]);
    si_djj_mean = mean(si_regrid(:,:,:,idx1),4);
    si_diff = si_djj_mean - si_obs_djj_mean;
    eval(['si_' model_abbs{i} '_djj_mean=si_djj_mean;']);
    eval(['si_' model_abbs{i} '_diff=si_diff;']);
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

% Plot surface silicate (Fig.5)
fig5=figure(1);
set(gcf,'Position',[100 100 1200 600])
sp0=subplot(3,6,1);
pos = get(gca, 'Position'); 
set(gca, 'Position', [pos(1) pos(2) pos(3)*1.15 pos(4)*1.15]); 
m_proj('stereographic','lat',-90,'long',0,'radius',60);
m_pcolor(lon2r,lat2r,si_obs_djj_mean(:,:,1))
shading flat;
colormap(sp0,flipud(m_colmap('blues')));
m_line(line_stf(:,1),line_stf(:,2),'color','k','LineStyle','--','linewi',.5)
m_line(line_saf(:,1),line_saf(:,2),'color','k','LineStyle','--','linewi',.5)
m_line(line_pf(:,1),line_pf(:,2),'color','k','LineStyle','--','linewi',.5)
m_coast('patch','w');
m_coast('linewidth',1,'color','k');
m_grid('xtick',6,'xticklabel',[],'tickdir','out','ytick',[],'linest','none','xaxisloc','top','yaxislocation','middle','backcolor',[0.7 0.7 0.7]);
title('WOA','FontSize',12)
clim([0 80])
c=colorbar;
c.Ticks = 0:20:80;
c.Label.String = 'Si (mmol/m^{3})';
c.Label.FontSize = 11;
set(c,'Location','southoutside','Position',[0.14, 0.65, 0.1, 0.02])

si_acc_diff = NaN(size(si_gfd_diff));
si_can_diff = NaN(size(si_gfd_diff));
si_mir_diff = NaN(size(si_gfd_diff));
plot_abbs = {'acc','can','ces','cmc','cnr','gfd','ips','mir','mhm','mhr','mlr','nlm','nmm','uke'};
plot_longname = {'ACCESS-ESM1-5^{*}','CanESM5^{*}','CESM2','CMCC-ESM2','CNRM-ESM2-1','GFDL-ESM4',...
    'IPSL-CM6A-LR','MIROC-ES2L^{*}','MPI-ESM-1-2-HAM','MPI-ESM1-2-HR',...
    'MPI-ESM1-2-LR','NorESM2-LM','NorESM2-MM','UKESM1-0-LL'};

for s = 1:length(splist)
    plotdata = eval(['si_' plot_abbs{s} '_diff']);
    sp=subplot(3,6,splist(s));
    pos = get(gca, 'Position'); 
    set(gca, 'Position', [pos(1) pos(2) pos(3)*1.15 pos(4)*1.15]); 
    m_proj('stereographic','lat',-90,'long',0,'radius',60);
    m_pcolor(lon2r,lat2r,plotdata(:,:,1))
    shading flat;
    colormap(sp,b2r(-80,80));
    m_line(line_stf(:,1),line_stf(:,2),'color','k','LineStyle','--','linewi',.5)
    m_line(line_saf(:,1),line_saf(:,2),'color','k','LineStyle','--','linewi',.5)
    m_line(line_pf(:,1),line_pf(:,2),'color','k','LineStyle','--','linewi',.5)
    m_coast('patch','w');
    m_coast('linewidth',1,'color','k');
    m_grid('xtick',6,'xticklabel',[],'tickdir','out','ytick',[],'linest','none','xaxisloc','top','yaxislocation','middle','backcolor',[0.7 0.7 0.7]);
    title(plot_longname{s},'FontSize',12)
end
c1=colorbar;
c1.Ticks = -80:40:80;
set(c1,'Location','eastoutside','Position',[0.82, 0.1, 0.01, 0.25])
c1.Label.String = 'Si (mmol/m^{3})';
c1.Label.FontSize = 11;

print(fig5,'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Figure\figs5','-djpeg','-r300');

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

mbe_si_os = nan(length(plot_abbs),5);
for i = 1:length(plot_abbs)
    si_diff_os = eval(['si_' plot_abbs{i} '_diff(:,:,1)']);
    mbe_si_os(i,1) = sum(si_diff_os(~isnan(si_diff_os)).*cosd(lat2r(~isnan(si_diff_os))))...
        ./sum((cosd(lat2r(~isnan(si_diff_os)))));
    mbe_si_os(i,2) = sum(si_diff_os(~isnan(si_diff_os)&idx_stz).*cosd(lat2r(~isnan(si_diff_os)&idx_stz)))...
        ./sum(cosd(lat2r(~isnan(si_diff_os)&idx_stz)));
    mbe_si_os(i,3) = sum(si_diff_os(~isnan(si_diff_os)&idx_saz).*cosd(lat2r(~isnan(si_diff_os)&idx_saz)))...
        ./sum(cosd(lat2r(~isnan(si_diff_os)&idx_saz)));
    mbe_si_os(i,4) = sum(si_diff_os(~isnan(si_diff_os)&idx_pfz).*cosd(lat2r(~isnan(si_diff_os)&idx_pfz)))...
        ./sum(cosd(lat2r(~isnan(si_diff_os)&idx_pfz)));
    mbe_si_os(i,5) = sum(si_diff_os(~isnan(si_diff_os)&idx_az).*cosd(lat2r(~isnan(si_diff_os)&idx_az)))...
        ./sum(cosd(lat2r(~isnan(si_diff_os)&idx_az)));
end

% Plot MBE (Fig.6)
col = [0.4 0 0.6; 1 0.7 0.7; 0 0.8 0; 1 0.8 0; 0.4 0.7 1];
fig6=figure(2);
set(gcf,'Position',[100,100,1200,500])
b=bar(mbe_si_os);
for c = 1:length(b)
    b(c).FaceColor = col(c,:);
end
ylabel('Mean Bias Error (mmol/m^{3})')
set(gca,'XTickLabel',plot_longname,'FontSize',12)
legend({'SO','STZ','SAZ','PFZ','AZ'})

print(fig6,'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Figure\fig6','-djpeg','-r300');

% Taylor Diagram
% Calculation of different parameters
tdlist = {'obs','acc','can','ces','cmc','cnr','gfd','ips','mir','mhm','mhr','mlr','nlm','nmm','uke'};
si_acc_djj_mean = NaN(size(si_obs_djj_mean));
si_can_djj_mean = NaN(size(si_obs_djj_mean));
si_mir_djj_mean = NaN(size(si_obs_djj_mean));
for i = 1:length(tdlist)
    si_so = eval(['si_' tdlist{i} '_djj_mean(:,:,1)']);
    si_stz = si_so(idx_stz);
    si_saz = si_so(idx_saz);
    si_pfz = si_so(idx_pfz);
    si_az = si_so(idx_az);
    eval(['si_so_' tdlist{i} '=si_so(:);']);
    eval(['si_stz_' tdlist{i} '=si_stz;']);
    eval(['si_saz_' tdlist{i} '=si_saz;']);
    eval(['si_pfz_' tdlist{i} '=si_pfz;']);
    eval(['si_az_' tdlist{i} '=si_az;']);
end

rglist = {'so','stz','saz','pfz','az'};

valid_idx_so = ~isnan(si_so_obs) & ~isnan(si_so_ces) & ...
    ~isnan(si_so_cmc) & ~isnan(si_so_gfd) & ~isnan(si_so_ips) & ...
    ~isnan(si_so_mhm) & ~isnan(si_so_mhr) & ~isnan(si_so_mlr) & ...
    ~isnan(si_so_nlm) & ~isnan(si_so_nmm);

valid_idx_stz = ~isnan(si_stz_obs) & ~isnan(si_stz_ces) & ...
    ~isnan(si_stz_cmc) & ~isnan(si_stz_gfd) & ~isnan(si_stz_ips) & ...
    ~isnan(si_stz_mhm) & ~isnan(si_stz_mhr) & ~isnan(si_stz_mlr) & ...
    ~isnan(si_stz_nlm) & ~isnan(si_stz_nmm);

valid_idx_saz = ~isnan(si_saz_obs) & ~isnan(si_saz_ces) & ...
    ~isnan(si_saz_cmc) & ~isnan(si_saz_gfd) & ~isnan(si_saz_ips) & ...
    ~isnan(si_saz_mhm) & ~isnan(si_saz_mhr) & ~isnan(si_saz_mlr) & ...
    ~isnan(si_saz_nlm) & ~isnan(si_saz_nmm);

valid_idx_pfz = ~isnan(si_pfz_obs) & ~isnan(si_pfz_ces) & ...
    ~isnan(si_pfz_cmc) & ~isnan(si_pfz_gfd) & ~isnan(si_pfz_ips) & ...
    ~isnan(si_pfz_mhm) & ~isnan(si_pfz_mhr) & ~isnan(si_pfz_mlr) & ...
    ~isnan(si_pfz_nlm) & ~isnan(si_pfz_nmm);

valid_idx_az = ~isnan(si_az_obs) & ~isnan(si_az_ces) & ...
    ~isnan(si_az_cmc) & ~isnan(si_az_gfd) & ~isnan(si_az_ips) & ...
    ~isnan(si_az_mhm) & ~isnan(si_az_mhr) & ~isnan(si_az_mlr) & ...
    ~isnan(si_az_nlm) & ~isnan(si_az_nmm);

for j = 1:length(rglist)
    for i = 1:length(tdlist)
        eval(['si_' rglist{j} '_' tdlist{i} '_nonan=si_' rglist{j} '_' tdlist{i} '(valid_idx_' rglist{j} ');']);
    end
end

stats = [];
for j = 1:length(rglist)
    for i = 1:length(plot_abbs)
        S = SStats(eval(['si_' rglist{j} '_obs_nonan']),eval(['si_' rglist{j} '_' plot_abbs{i} '_nonan']));
        stats(j,:,i+1) = S;
    end
    stats(j,:,1) = SStats(eval(['si_' rglist{j} '_obs_nonan']),eval(['si_' rglist{j} '_obs_nonan']));

    stdev = stats(j,2,1); 

    stats(j,2,:) = stats(j,2,:)/stdev; 
    stats(j,3,:) = stats(j,3,:)/stdev;
end

stats(:,1,:) = [];

% Taylor diagram for suraface silicate (Fig.S3)
td_longname = {' ','ACCESS-ESM1-5^{*}','CanESM5^{*}','CESM2','CMCC-ESM2','CNRM-ESM2-1','GFDL-ESM4',...
    'IPSL-CM6A-LR','MIROC-ES2L^{*}','MPI-ESM-1-2-HAM','MPI-ESM1-2-HR',...
    'MPI-ESM1-2-LR','NorESM2-LM','NorESM2-MM','UKESM1-0-LL'};

figs3=figure(3);
set(gcf,'Position',[100 100 1500 1000])
subplot(2,3,1)
[hp, ht, axl] = taylor_diagram(squeeze(stats(1,1,:)),squeeze(stats(1,2,:)),squeeze(stats(1,3,:)), ...
    'markerLabel',td_longname,'markerLegend', 'on', ...
    'styleSTD', '-', 'colOBS','k', 'markerObs','o', ...
    'markerSize',6, 'tickRMS',0:0.5:1, ...
    'tickRMSangle', 115, 'showlabelsRMS', 'on', ...
    'titleRMS','on', 'titleOBS','Obs');
t1=title('Surface silicate in SO','FontSize',14);
t1.Position(2) = t1.Position(2) + 0.12;

subplot(2,3,2)
[hp, ht, axl] = taylor_diagram(squeeze(stats(2,1,:)),squeeze(stats(2,2,:)),squeeze(stats(2,3,:)), ...
    'markerLabel',td_longname, 'markerLegend', 'on', ...
    'styleSTD', '-', 'colOBS','k', 'markerObs','o', ...
    'markerSize',6, 'tickRMS',0:3:15, ...
    'tickRMSangle', 115, 'showlabelsRMS', 'on', ...
    'titleRMS','on', 'titleOBS','Obs');
t2=title('Surface silicate in STZ','FontSize',14);
t2.Position(2) = t2.Position(2) + 10;

subplot(2,3,3)
[hp, ht, axl] = taylor_diagram(squeeze(stats(3,1,:)),squeeze(stats(3,2,:)),squeeze(stats(3,3,:)), ...
    'markerLabel',td_longname, 'markerLegend', 'on', ...
    'styleSTD', '-', 'colOBS','k', 'markerObs','o', ...
    'markerSize',6, 'tickRMS',0:2:10, ...
    'tickRMSangle', 115, 'showlabelsRMS', 'on', ...
    'titleRMS','on', 'titleOBS','Obs');
t3=title('Surface silicate in SAZ','FontSize',14);
t3.Position(2) = t3.Position(2) + 0.5;

subplot(2,3,4)
[hp, ht, axl] = taylor_diagram(squeeze(stats(4,1,:)),squeeze(stats(4,2,:)),squeeze(stats(4,3,:)), ...
    'markerLabel',td_longname, 'markerLegend', 'on', ...
    'styleSTD', '-', 'colOBS','k', 'markerObs','o', ...
    'markerSize',6, 'tickRMS',0:1:2, ...
    'tickRMSangle', 115, 'showlabelsRMS', 'on', ...
    'titleRMS','on', 'titleOBS','Obs');
t4=title('Surface silicate in PFZ','FontSize',14);
t4.Position(2) = t4.Position(2) + 0.2;

subplot(2,3,5)
[hp, ht, axl] = taylor_diagram(squeeze(stats(5,1,:)),squeeze(stats(5,2,:)),squeeze(stats(5,3,:)), ...
    'markerLabel',td_longname, 'markerLegend', 'on', ...
    'styleSTD', '-', 'colOBS','k', 'markerObs','o', ...
    'markerSize',6, 'tickRMS',0:0.5:1, ...
    'tickRMSangle', 115, 'showlabelsRMS', 'on', ...
    'titleRMS','on', 'titleOBS','Obs');
t5=title('Surface silicate in AZ','FontSize',14);
t5.Position(2) = t5.Position(2) + 0.12;

print(figs3,'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Figure\figs3','-djpeg','-r300');
