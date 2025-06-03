% Evaluating the performance of CMIP6 models in simulating Southern Ocean
% biogeochemistry
% Author: Cheng, M., Ellwood, M., and Maher, N. 
% Last modified: 26/05/2025
% POC related analysis (Section 3.3)
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

pocfile_obs = 'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Regrided_data\Observation\plankton_copernicus_final.nc';

phycfile_acc = 'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Regrided_data\CMIP\phyc_ACCESS-ESM1-5_final.nc';
phycfile_can = 'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Regrided_data\CMIP\phyc_CanESM5_final.nc';
phycfile_ces = 'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Regrided_data\CMIP\phyc_CESM2_final.nc';
phycfile_cmc = 'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Regrided_data\CMIP\phyc_CMCC-ESM2_final.nc';
phycfile_cnr = 'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Regrided_data\CMIP\phyc_CNRM-ESM2-1_final.nc';
phycfile_gfd = 'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Regrided_data\CMIP\phyc_GFDL-ESM4_final.nc';
phycfile_ips = 'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Regrided_data\CMIP\phyc_IPSL-CM6A-LR_final.nc';
phycfile_mir = 'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Regrided_data\CMIP\phyc_MIROC-ES2L_final.nc';
phycfile_mhm = 'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Regrided_data\CMIP\phyc_MPI-ESM-1-2-HAM_final.nc';
phycfile_mhr = 'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Regrided_data\CMIP\phyc_MPI-ESM1-2-HR_final.nc';
phycfile_mlr = 'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Regrided_data\CMIP\phyc_MPI-ESM1-2-LR_final.nc';
phycfile_nlm = 'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Regrided_data\CMIP\phyc_NorESM2-LM_final.nc';
phycfile_nmm = 'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Regrided_data\CMIP\phyc_NorESM2-MM_final.nc';
phycfile_uke = 'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Regrided_data\CMIP\phyc_UKESM1-0-LL_final.nc';

detocfile_acc = 'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Regrided_data\CMIP\detoc_ACCESS-ESM1-5_final.nc';
detocfile_can = 'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Regrided_data\CMIP\detoc_CanESM5_final.nc';
detocfile_cmc = 'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Regrided_data\CMIP\detoc_CMCC-ESM2_final.nc';
detocfile_cnr = 'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Regrided_data\CMIP\detoc_CNRM-ESM2-1_final.nc';
detocfile_gfd = 'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Regrided_data\CMIP\detoc_GFDL-ESM4_final.nc';
detocfile_ips = 'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Regrided_data\CMIP\detoc_IPSL-CM6A-LR_final.nc';
detocfile_mhm = 'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Regrided_data\CMIP\detoc_MPI-ESM-1-2-HAM_final.nc';
detocfile_mhr = 'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Regrided_data\CMIP\detoc_MPI-ESM1-2-HR_final.nc';
detocfile_mlr = 'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Regrided_data\CMIP\detoc_MPI-ESM1-2-LR_final.nc';
detocfile_nlm = 'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Regrided_data\CMIP\detoc_NorESM2-LM_final.nc';
detocfile_nmm = 'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Regrided_data\CMIP\detoc_NorESM2-MM_final.nc';
detocfile_uke = 'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Regrided_data\CMIP\detoc_UKESM1-0-LL_final.nc';

zoocfile_acc = 'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Regrided_data\CMIP\zooc_ACCESS-ESM1-5_final.nc';
zoocfile_can = 'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Regrided_data\CMIP\zooc_CanESM5_final.nc';
zoocfile_ces = 'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Regrided_data\CMIP\zooc_CESM2_final.nc';
zoocfile_cmc = 'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Regrided_data\CMIP\zooc_CMCC-ESM2_final.nc';
zoocfile_cnr = 'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Regrided_data\CMIP\zooc_CNRM-ESM2-1_final.nc';
zoocfile_gfd = 'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Regrided_data\CMIP\zooc_GFDL-ESM4_final.nc';
zoocfile_ips = 'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Regrided_data\CMIP\zooc_IPSL-CM6A-LR_final.nc';
zoocfile_mir = 'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Regrided_data\CMIP\zooc_MIROC-ES2L_final.nc';
zoocfile_mhm = 'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Regrided_data\CMIP\zooc_MPI-ESM-1-2-HAM_final.nc';
zoocfile_mhr = 'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Regrided_data\CMIP\zooc_MPI-ESM1-2-HR_final.nc';
zoocfile_mlr = 'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Regrided_data\CMIP\zooc_MPI-ESM1-2-LR_final.nc';
zoocfile_nlm = 'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Regrided_data\CMIP\zooc_NorESM2-LM_final.nc';
zoocfile_nmm = 'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Regrided_data\CMIP\zooc_NorESM2-MM_final.nc';
zoocfile_uke = 'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Regrided_data\CMIP\zooc_UKESM1-0-LL_final.nc';

baccfile_cmc = 'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Regrided_data\CMIP\bacc_CMCC-ESM2_final.nc';
baccfile_gfd = 'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Regrided_data\CMIP\bacc_GFDL-ESM4_final.nc';

model_abbs = {'acc','can','ces','cmc','cnr','gfd','ips','mir','mhm','mhr','mlr','nlm','nmm','uke'};
model_longname = {'ACCESS-ESM1-5','CanESM5','CESM2','CMCC-ESM2','CNRM-ESM2-1','GFDL-ESM4',...
    'IPSL-CM6A-LR','MIROC-ES2L','MPI-ESM-1-2-HAM','MPI-ESM1-2-HR',...
    'MPI-ESM1-2-LR','NorESM2-LM','NorESM2-MM','UKESM1-0-LL'};
chl_lev = {'lev','lev','lev_partial','lev','lev','lev','olevel','lev','lev','lev','lev','lev','lev','lev'};

all_abbs = {'obs','acc','can','ces','cmc','cnr','gfd','ips','mir','mhm','mhr','mlr','nlm','nmm','uke'};
all_longname = {'Copernicus','ACCESS-ESM1-5','CanESM5','CESM2','CMCC-ESM2','CNRM-ESM2-1','GFDL-ESM4',...
    'IPSL-CM6A-LR','MIROC-ES2L','MPI-ESM-1-2-HAM','MPI-ESM1-2-HR',...
    'MPI-ESM1-2-LR','NorESM2-LM','NorESM2-MM','UKESM1-0-LL'};

poc_lev = {'st_ocean','lev','lev','lev','lev','lev','olevel','lev','lev','lev','lev','lev','lev','lev'};

detoc_abbs = {'acc','can','cmc','cnr','gfd','ips','mhm','mhr','mlr','nlm','nmm','uke'};
detoc_lev = {'st_ocean','lev','lev','lev','lev','olevel','lev','lev','lev','lev','lev','lev'};

bac_abbs = {'cmc','gfd'};

lon = ncread(pocfile_obs,'lon');
lat = ncread(pocfile_obs,'lat');
lev = ncread(pocfile_obs,'depth');
tim = ncread(pocfile_obs,'time');
chl_obs_read = ncread(chlfile_obs,'chl');
poc_obs_read = ncread(pocfile_obs,'poc');

lonr = [lon;lon(1)];

[lat2,lon2] = meshgrid(lat,lon);
[lat2r,lon2r] = meshgrid(lat,lonr);

chl_obs = chl_obs_read([1:360,1],:,2:end,:);
poc_obs = poc_obs_read([1:360,1],:,2:end,:);

% Calculate mean in austral summer (DJJ) from 2000 to 2014
idx = ([1,2,12].*ones(15,1)+(0:12:14*12)')';
idx1 = idx(:);
chl_obs_djj_mean = nanmean(chl_obs(:,:,:,idx1),4);
poc_obs_djj_mean = nanmean(poc_obs,4);

for i = 1:length(model_abbs)
    filename = eval(['chlfile_' model_abbs{i}]);
    chl_read = ncread(filename,'chl');
    lev_read = ncread(filename,chl_lev{i});
    eval(['chl_' model_abbs{i} '_read=chl_read;']);
    eval(['lev_' model_abbs{i} '=lev_read;']);
end

% Optimise some depth mannually
lev = double(lev(2:end));
lev_mhm(1) = 5;
lev_mhr(1) = 5;
lev_mlr(1) = 5;

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
    clear(['chl_' model_abbs{i} '_read']);
end

chl_cnr = chl_cnr*1e-3;
chl_ips = chl_ips*1e-3;

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

for i = 1:length(model_abbs)
    chl_regrid = eval(['chl_' model_abbs{i}]);
    chl_djj_mean = mean(chl_regrid(:,:,:,idx1),4);
    chl_diff = chl_djj_mean - chl_obs_djj_mean;
    eval(['chl_' model_abbs{i} '_djj_mean=chl_djj_mean;']);
    eval(['chl_' model_abbs{i} '_diff=chl_diff;']);
    clear(['chl_' model_abbs{i}]);
end

for i = 1:length(model_abbs)
    filename = eval(['phycfile_' model_abbs{i}]);
    phyc_read = ncread(filename,'phyc');
    lev_read = ncread(filename,poc_lev{i});
    eval(['phyc_' model_abbs{i} '_read=phyc_read;']);
    eval(['lev_' model_abbs{i} '=lev_read;']);
end

lev_ces = lev_ces/100;
lev_mhm(1) = 5;
lev_mhr(1) = 5;
lev_mlr(1) = 5;

for i = 1:length(model_abbs)
    phyc_regrid = NaN(size(chl_obs_read,1),size(chl_obs_read,2),size(chl_obs_read,3)-1,15);
    phyc_read = eval(['phyc_' model_abbs{i} '_read']);
    lev_read = eval(['lev_' model_abbs{i}]);
    for ji = 1:size(phyc_regrid,1)
        for jj = 1:size(phyc_regrid,2)
            for jk = 1:size(phyc_regrid,4)
                phyc_regrid(ji,jj,:,jk) = interp1(lev_read,squeeze(phyc_read(ji,jj,:,jk)),lev)*12000;
            end
        end
    end
    phyc361 = phyc_regrid([1:360,1],:,:,:);
    eval(['phyc_' model_abbs{i} '=phyc361;']);
    clear(['phyc_' model_abbs{i} '_read']);
end

phyc_acc = phyc_acc*6.625e-3;

phyc_can_regrid = phyc_can;
phyc_cmc_regrid = phyc_cmc;
phyc_cnr_regrid = phyc_cnr;
phyc_ips_regrid = phyc_ips;
phyc_mhm_regrid = phyc_mhm;
phyc_mlr_regrid = phyc_mlr;

phyc_can(74,:,:,:) = mean(phyc_can_regrid([73,75],:,:,:));
phyc_cmc(73,:,:,:) = mean(phyc_cmc_regrid([72,75],:,:,:));
phyc_cmc(74,:,:,:) = mean(phyc_cmc_regrid([72,75],:,:,:));
phyc_cnr(73,:,:,:) = mean(phyc_cnr_regrid([72,75],:,:,:));
phyc_cnr(74,:,:,:) = mean(phyc_cnr_regrid([72,75],:,:,:));
phyc_ips(73,:,:,:) = mean(phyc_ips_regrid([72,75],:,:,:));
phyc_ips(74,:,:,:) = mean(phyc_ips_regrid([72,75],:,:,:));

lon_fill_mhm = [147;148;148;149;149;150;150;151;151;152;152;153;153;153;...
    154;154;154;155;155;155;155;156;156;156;156;157;157;157;157;157];
lat_fill_mhm = [24;25;26;27;28;29;30;31;32;34;35;37;38;39;41;42;43;45;46;...
    47;48;50;51;52;53;56;57;58;59;60];

for f = 1:length(lon_fill_mhm)
    for ji = 1:length(lon_fill_mhm)
        for jj = 1:size(phyc_mhm,3)
            for jk = 1:size(phyc_mhm,4)
                phyc_mhm(lon_fill_mhm(ji),lat_fill_mhm(ji),jj,jk) = ...
                    mean(phyc_mhm_regrid([lon_fill_mhm(ji)-1,lon_fill_mhm(ji)+1],lat_fill_mhm(ji),jj,jk));
            end
        end
    end
end

for f = 1:length(lon_fill_mhm)
    for ji = 1:length(lon_fill_mhm)
        for jj = 1:size(phyc_mlr,3)
            for jk = 1:size(phyc_mlr,4)
                phyc_mlr(lon_fill_mhm(ji),lat_fill_mhm(ji),jj,jk) = ...
                    mean(phyc_mlr_regrid([lon_fill_mhm(ji)-1,lon_fill_mhm(ji)+1],lat_fill_mhm(ji),jj,jk));
            end
        end
    end
end

for i = 1:length(model_abbs)
    phyc_regrid = eval(['phyc_' model_abbs{i}]);
    phyc_djj_mean = mean(phyc_regrid,4);
    eval(['phyc_' model_abbs{i} '_djj_mean=phyc_djj_mean;']);
    clear(['phyc_' model_abbs{i}]);
end

for i = 1:length(detoc_abbs)
    filename = eval(['detocfile_' detoc_abbs{i}]);
    detoc_read = ncread(filename,'detoc');
    lev_read = ncread(filename,detoc_lev{i});
    eval(['detoc_' detoc_abbs{i} '_read=detoc_read;']);
    eval(['lev_' detoc_abbs{i} '=lev_read;']);
end

lev_mhm(1) = 5;
lev_mhr(1) = 5;
lev_mlr(1) = 5;

for i = 1:length(detoc_abbs)
    detoc_regrid = NaN(size(chl_obs_read,1),size(chl_obs_read,2),size(chl_obs_read,3)-1,15);
    detoc_read = eval(['detoc_' detoc_abbs{i} '_read']);
    lev_read = eval(['lev_' detoc_abbs{i}]);
    for ji = 1:size(detoc_regrid,1)
        for jj = 1:size(detoc_regrid,2)
            for jk = 1:size(detoc_regrid,4)
                detoc_regrid(ji,jj,:,jk) = interp1(lev_read,squeeze(detoc_read(ji,jj,:,jk)),lev)*12000;
            end
        end
    end
    detoc361 = detoc_regrid([1:360,1],:,:,:);
    eval(['detoc_' detoc_abbs{i} '=detoc361;']);
    clear(['detoc_' detoc_abbs{i} '_read']);
end

detoc_acc = detoc_acc*6.625e-3;

detoc_can_regrid = detoc_can;
detoc_cmc_regrid = detoc_cmc;
detoc_cnr_regrid = detoc_cnr;
detoc_ips_regrid = detoc_ips;
detoc_mhm_regrid = detoc_mhm;
detoc_mlr_regrid = detoc_mlr;

detoc_can(74,:,:,:) = mean(detoc_can_regrid([73,75],:,:,:));
detoc_cmc(73,:,:,:) = mean(detoc_cmc_regrid([72,75],:,:,:));
detoc_cmc(74,:,:,:) = mean(detoc_cmc_regrid([72,75],:,:,:));
detoc_cnr(73,:,:,:) = mean(detoc_cnr_regrid([72,75],:,:,:));
detoc_cnr(74,:,:,:) = mean(detoc_cnr_regrid([72,75],:,:,:));
detoc_ips(73,:,:,:) = mean(detoc_ips_regrid([72,75],:,:,:));
detoc_ips(74,:,:,:) = mean(detoc_ips_regrid([72,75],:,:,:));

lon_fill_mhm = [147;148;148;149;149;150;150;151;151;152;152;153;153;153;...
    154;154;154;155;155;155;155;156;156;156;156;157;157;157;157;157];
lat_fill_mhm = [24;25;26;27;28;29;30;31;32;34;35;37;38;39;41;42;43;45;46;...
    47;48;50;51;52;53;56;57;58;59;60];

for f = 1:length(lon_fill_mhm)
    for ji = 1:length(lon_fill_mhm)
        for jj = 1:size(detoc_mhm,3)
            for jk = 1:size(detoc_mhm,4)
                detoc_mhm(lon_fill_mhm(ji),lat_fill_mhm(ji),jj,jk) = ...
                    mean(detoc_mhm_regrid([lon_fill_mhm(ji)-1,lon_fill_mhm(ji)+1],lat_fill_mhm(ji),jj,jk));
            end
        end
    end
end

for f = 1:length(lon_fill_mhm)
    for ji = 1:length(lon_fill_mhm)
        for jj = 1:size(detoc_mlr,3)
            for jk = 1:size(detoc_mlr,4)
                detoc_mlr(lon_fill_mhm(ji),lat_fill_mhm(ji),jj,jk) = ...
                    mean(detoc_mlr_regrid([lon_fill_mhm(ji)-1,lon_fill_mhm(ji)+1],lat_fill_mhm(ji),jj,jk));
            end
        end
    end
end

for i = 1:length(detoc_abbs)
    detoc_regrid = eval(['detoc_' detoc_abbs{i}]);
    detoc_djj_mean = mean(detoc_regrid,4);
    eval(['detoc_' detoc_abbs{i} '_djj_mean=detoc_djj_mean;']);
    clear(['detoc_' detoc_abbs{i}]);
end

detoc_ces_djj_mean = zeros(size(detoc_djj_mean));
detoc_mir_djj_mean = zeros(size(detoc_djj_mean));

for i = 1:length(model_abbs)
    filename = eval(['zoocfile_' model_abbs{i}]);
    zooc_read = ncread(filename,'zooc');
    lev_read = ncread(filename,poc_lev{i});
    eval(['zooc_' model_abbs{i} '_read=zooc_read;']);
    eval(['lev_' model_abbs{i} '=lev_read;']);
end

lev_ces = lev_ces/100;
lev_mhm(1) = 5;
lev_mhr(1) = 5;
lev_mlr(1) = 5;

for i = 1:length(model_abbs)
    zooc_regrid = NaN(size(chl_obs_read,1),size(chl_obs_read,2),size(chl_obs_read,3)-1,15);
    zooc_read = eval(['zooc_' model_abbs{i} '_read']);
    lev_read = eval(['lev_' model_abbs{i}]);
    for ji = 1:size(zooc_regrid,1)
        for jj = 1:size(zooc_regrid,2)
            for jk = 1:size(zooc_regrid,4)
                zooc_regrid(ji,jj,:,jk) = interp1(lev_read,squeeze(zooc_read(ji,jj,:,jk)),lev)*12000;
            end
        end
    end
    zooc361 = zooc_regrid([1:360,1],:,:,:);
    eval(['zooc_' model_abbs{i} '=zooc361;']);
    clear(['zooc_' model_abbs{i} '_read']);
end

zooc_acc = zooc_acc*6.625e-3;

zooc_can_regrid = zooc_can;
zooc_cmc_regrid = zooc_cmc;
zooc_cnr_regrid = zooc_cnr;
zooc_ips_regrid = zooc_ips;
zooc_mhm_regrid = zooc_mhm;
zooc_mlr_regrid = zooc_mlr;

zooc_can(74,:,:,:) = mean(zooc_can_regrid([73,75],:,:,:));
zooc_cmc(73,:,:,:) = mean(zooc_cmc_regrid([72,75],:,:,:));
zooc_cmc(74,:,:,:) = mean(zooc_cmc_regrid([72,75],:,:,:));
zooc_cnr(73,:,:,:) = mean(zooc_cnr_regrid([72,75],:,:,:));
zooc_cnr(74,:,:,:) = mean(zooc_cnr_regrid([72,75],:,:,:));
zooc_ips(73,:,:,:) = mean(zooc_ips_regrid([72,75],:,:,:));
zooc_ips(74,:,:,:) = mean(zooc_ips_regrid([72,75],:,:,:));

lon_fill_mhm = [147;148;148;149;149;150;150;151;151;152;152;153;153;153;...
    154;154;154;155;155;155;155;156;156;156;156;157;157;157;157;157];
lat_fill_mhm = [24;25;26;27;28;29;30;31;32;34;35;37;38;39;41;42;43;45;46;...
    47;48;50;51;52;53;56;57;58;59;60];

for f = 1:length(lon_fill_mhm)
    for ji = 1:length(lon_fill_mhm)
        for jj = 1:size(zooc_mhm,3)
            for jk = 1:size(zooc_mhm,4)
                zooc_mhm(lon_fill_mhm(ji),lat_fill_mhm(ji),jj,jk) = ...
                    mean(zooc_mhm_regrid([lon_fill_mhm(ji)-1,lon_fill_mhm(ji)+1],lat_fill_mhm(ji),jj,jk));
            end
        end
    end
end

for f = 1:length(lon_fill_mhm)
    for ji = 1:length(lon_fill_mhm)
        for jj = 1:size(zooc_mlr,3)
            for jk = 1:size(zooc_mlr,4)
                zooc_mlr(lon_fill_mhm(ji),lat_fill_mhm(ji),jj,jk) = ...
                    mean(zooc_mlr_regrid([lon_fill_mhm(ji)-1,lon_fill_mhm(ji)+1],lat_fill_mhm(ji),jj,jk));
            end
        end
    end
end

for i = 1:length(model_abbs)
    zooc_regrid = eval(['zooc_' model_abbs{i}]);
    zooc_djj_mean = mean(zooc_regrid,4);
    eval(['zooc_' model_abbs{i} '_djj_mean=zooc_djj_mean;']);
    clear(['zooc_' model_abbs{i}]);
end

for i = 1:length(bac_abbs)
    filename = eval(['baccfile_' bac_abbs{i}]);
    bacc_read = ncread(filename,'bacc');
    lev_read = ncread(filename,'lev');
    eval(['bacc_' bac_abbs{i} '_read=bacc_read;']);
    eval(['lev_' bac_abbs{i} '=lev_read;']);
end

for i = 1:length(bac_abbs)
    bacc_regrid = NaN(size(chl_obs_read,1),size(chl_obs_read,2),size(chl_obs_read,3)-1,15);
    bacc_read = eval(['bacc_' bac_abbs{i} '_read']);
    lev_read = eval(['lev_' bac_abbs{i}]);
    for ji = 1:size(bacc_regrid,1)
        for jj = 1:size(bacc_regrid,2)
            for jk = 1:size(bacc_regrid,4)
                bacc_regrid(ji,jj,:,jk) = interp1(lev_read,squeeze(bacc_read(ji,jj,:,jk)),lev)*12000;
            end
        end
    end
    bacc361 = bacc_regrid([1:360,1],:,:,:);
    eval(['bacc_' bac_abbs{i} '=bacc361;']);
    clear(['bacc_' bac_abbs{i} '_read']);
end

bacc_cmc_regrid = bacc_cmc;

bacc_cmc(73,:,:,:) = mean(bacc_cmc_regrid([72,75],:,:,:));
bacc_cmc(74,:,:,:) = mean(bacc_cmc_regrid([72,75],:,:,:));

for i = 1:length(bac_abbs)
    bacc_regrid = eval(['bacc_' bac_abbs{i}]);
    bacc_djj_mean = mean(bacc_regrid,4);
    eval(['bacc_' bac_abbs{i} '_djj_mean=bacc_djj_mean;']);
    clear(['bacc_' bac_abbs{i}]);
end

for i = 1:length(model_abbs)
    if ismember(model_abbs{i},bac_abbs)
        eval(['poc_' model_abbs{i} '_djj_mean=phyc_' model_abbs{i} '_djj_mean+detoc_' model_abbs{i} '_djj_mean+zooc_' model_abbs{i} '_djj_mean+bacc_' model_abbs{i} '_djj_mean;']);
    else
        eval(['poc_' model_abbs{i} '_djj_mean=phyc_' model_abbs{i} '_djj_mean+detoc_' model_abbs{i} '_djj_mean+zooc_' model_abbs{i} '_djj_mean;']);
    end
    eval(['poc_' model_abbs{i} '_diff=poc_' model_abbs{i} '_djj_mean-poc_obs_djj_mean;']);
end

% Calculate the vertical integration
lev2 = [0;lev];
chlz_obs = NaN(size(chl_obs_djj_mean,1),size(chl_obs_djj_mean,2));
chl_obs_djj_mean2 = NaN(size(chl_obs_djj_mean,1),size(chl_obs_djj_mean,2),length(lev2));
chl_obs_djj_mean2(:,:,1) = chl_obs_djj_mean(:,:,1);
chl_obs_djj_mean2(:,:,2:end) = chl_obs_djj_mean;
chl_obs_djj_mean2(isnan(chl_obs_djj_mean2)) = 0;
for ji = 1:size(chlz_obs,1)
    for jj = 1:size(chlz_obs,2)
        chlz_obs(ji,jj) = trapz(lev2(1:18),chl_obs_djj_mean2(ji,jj,1:18));
    end
end

chlz_obsp = chlz_obs;
chlz_obsp(chlz_obsp==0) = NaN;

for i = 1:length(model_abbs)
    chl_djj_mean = eval(['chl_' model_abbs{i} '_djj_mean']);
    chl_djj_mean2 = NaN(size(chl_djj_mean,1),size(chl_djj_mean,2),length(lev2));
    chl_djj_mean2(:,:,1) = chl_djj_mean(:,:,1);
    chl_djj_mean2(:,:,2:end) = chl_djj_mean;
    chl_djj_mean2(isnan(chl_djj_mean2)) = 0;
    chlz = NaN(size(chl_djj_mean,1),size(chl_djj_mean,2));
    for ji = 1:size(chl_djj_mean,1)
        for jj = 1:size(chl_djj_mean,2)
            chlz(ji,jj) = trapz(lev2(1:18),chl_djj_mean2(ji,jj,1:18));
        end
    end
    eval(['chlz_' model_abbs{i} '=chlz;']);
    eval(['chlz_' model_abbs{i} '_diff=chlz-chlz_obsp;'])
end

pocz_obs = NaN(size(poc_obs_djj_mean,1),size(poc_obs_djj_mean,2));
poc_obs_djj_mean2 = NaN(size(poc_obs_djj_mean,1),size(poc_obs_djj_mean,2),length(lev2));
poc_obs_djj_mean2(:,:,1) = poc_obs_djj_mean(:,:,1);
poc_obs_djj_mean2(:,:,2:end) = poc_obs_djj_mean;
poc_obs_djj_mean2(isnan(poc_obs_djj_mean2)) = 0;
for ji = 1:size(pocz_obs,1)
    for jj = 1:size(pocz_obs,2)
        pocz_obs(ji,jj) = trapz(lev2(1:18),poc_obs_djj_mean2(ji,jj,1:18));
    end
end

pocz_obsp = pocz_obs;
pocz_obsp(pocz_obsp==0) = NaN;

for i = 1:length(model_abbs)
    poc_djj_mean = eval(['poc_' model_abbs{i} '_djj_mean']);
    poc_djj_mean2 = NaN(size(poc_djj_mean,1),size(poc_djj_mean,2),length(lev2));
    poc_djj_mean2(:,:,1) = poc_djj_mean(:,:,1);
    poc_djj_mean2(:,:,2:end) = poc_djj_mean;
    poc_djj_mean2(isnan(poc_djj_mean2)) = 0;
    pocz = NaN(size(poc_djj_mean,1),size(poc_djj_mean,2));
    for ji = 1:size(poc_djj_mean,1)
        for jj = 1:size(poc_djj_mean,2)
            pocz(ji,jj) = trapz(lev2(1:18),poc_djj_mean2(ji,jj,1:18));
        end
    end
    pocz(isnan(pocz)) = 0;
    eval(['pocz_' model_abbs{i} '=pocz;']);
    eval(['pocz_' model_abbs{i} '_diff=pocz-pocz_obsp;'])
end

% Read front data
frontname = {'stf','saf','pf'};
for i = 1:length(frontname)
    txtfile = ['D:\OneDrive - Australian National University\PhD\Project1\MidTerm\Data\' frontname{i} '.txt'];
    % txtfile = ['C:\Users\u7392727\OneDrive - Australian National University\PhD\Project1\MidTerm\Data\' frontname{i} '.txt'];
    fid = fopen(txtfile,'r');
    header = fgetl(fid);
    data = textscan(fid,'%f %f');
    fclose(fid);
    line=[];
    line(:,1) = data{1};
    line(:,2) = data{2};
    eval(['line_' frontname{i} '=line;']);
end

splist = [2;3;4;5;6;8;9;10;11;12;14;15;16;17];
% Plot surface POC (Fig.11)
fig11=figure(1);
set(gcf,'Position',[100 100 1200 600])
sp0=subplot(3,6,1);
pos = get(gca, 'Position'); 
set(gca, 'Position', [pos(1) pos(2) pos(3)*1.15 pos(4)*1.15]); 
m_proj('stereographic','lat',-90,'long',0,'radius',60);
m_pcolor(lon2r,lat2r,poc_obs_djj_mean(:,:,1))
shading flat;
colormap(sp0,flipud(m_colmap('blue')));
m_line(line_stf(:,1),line_stf(:,2),'color','k','LineStyle','--','linewi',.5)
m_line(line_saf(:,1),line_saf(:,2),'color','k','LineStyle','--','linewi',.5)
m_line(line_pf(:,1),line_pf(:,2),'color','k','LineStyle','--','linewi',.5)
m_coast('patch','w');
m_coast('linewidth',1,'color','k');
m_grid('xtick',6,'xticklabel',[],'tickdir','out','ytick',[],'linest','none','xaxisloc','top','yaxislocation','middle','backcolor',[0.7 0.7 0.7]);
title('Copernicus','FontSize',12)
clim([0 150])
c=colorbar;
c.Label.String = 'POC (mg/m^{3})';
c.Label.FontSize = 11;
set(c,'Location','southoutside','Position',[0.14, 0.65, 0.1, 0.02])

for s = 1:length(model_abbs)
    plotdata = eval(['poc_' model_abbs{s} '_diff(:,:,1)']);
    sp=subplot(3,6,splist(s));
    pos = get(gca, 'Position'); 
    set(gca, 'Position', [pos(1) pos(2) pos(3)*1.15 pos(4)*1.15]); 
    m_proj('stereographic','lat',-90,'long',0,'radius',60);
    m_pcolor(lon2r,lat2r,plotdata(:,:,1))
    shading flat;
    colormap(sp,b2r(-150,150));
    m_line(line_stf(:,1),line_stf(:,2),'color','k','LineStyle','--','linewi',.5)
    m_line(line_saf(:,1),line_saf(:,2),'color','k','LineStyle','--','linewi',.5)
    m_line(line_pf(:,1),line_pf(:,2),'color','k','LineStyle','--','linewi',.5)
    m_coast('patch','w');
    m_coast('linewidth',1,'color','k');
    m_grid('xtick',6,'xticklabel',[],'tickdir','out','ytick',[],'linest','none','xaxisloc','top','yaxislocation','middle','backcolor',[0.7 0.7 0.7]);
    title(model_longname{s},'FontSize',12)
end
c1=colorbar;
set(c1,'Location','eastoutside','Position',[0.82, 0.1, 0.01, 0.25])
c1.Label.String = 'POC (mg/m^{3})';
c1.Label.FontSize = 11;

print(fig11,'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Figure\fig11','-djpeg','-r300');

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

mbe_poc = nan(length(model_abbs),5);
for i = 1:length(model_abbs)
    poc_diff = eval(['poc_' model_abbs{i} '_diff(:,:,1)']);
    mbe_poc(i,1) = sum(poc_diff(~isnan(poc_diff)).*cosd(lat2r(~isnan(poc_diff))))...
        ./sum((cosd(lat2r(~isnan(poc_diff)))));
    mbe_poc(i,2) = sum(poc_diff(~isnan(poc_diff)&idx_stz).*cosd(lat2r(~isnan(poc_diff)&idx_stz)))...
        ./sum(cosd(lat2r(~isnan(poc_diff)&idx_stz)));
    mbe_poc(i,3) = sum(poc_diff(~isnan(poc_diff)&idx_saz).*cosd(lat2r(~isnan(poc_diff)&idx_saz)))...
        ./sum(cosd(lat2r(~isnan(poc_diff)&idx_saz)));
    mbe_poc(i,4) = sum(poc_diff(~isnan(poc_diff)&idx_pfz).*cosd(lat2r(~isnan(poc_diff)&idx_pfz)))...
        ./sum(cosd(lat2r(~isnan(poc_diff)&idx_pfz)));
    mbe_poc(i,5) = sum(poc_diff(~isnan(poc_diff)&idx_az).*cosd(lat2r(~isnan(poc_diff)&idx_az)))...
        ./sum(cosd(lat2r(~isnan(poc_diff)&idx_az)));
end

% Plot MBE for surface POC (Fig.12)
col = [0.4 0 0.6; 1 0.7 0.7; 0 0.8 0; 1 0.8 0; 0.4 0.7 1];
fig12=figure(2);
set(gcf,'Position',[100,100,1200,500])
b=bar(mbe_poc);
for c = 1:length(b)
    b(c).FaceColor = col(c,:);
end
ylabel('Mean Bias Error (mg/m^{3})')
set(gca,'XTickLabel',model_longname,'FontSize',12)
legend({'SO','STZ','SAZ','PFZ','AZ'},'Location','northeast')

print(fig12,'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Figure\fig12','-djpeg','-r300');

% Taylor Diagram
% Calculation of different parameters
for i = 1:length(all_abbs)
    poc_so = eval(['poc_' all_abbs{i} '_djj_mean(:,:,1)']);
    poc_stz = poc_so(idx_stz);
    poc_saz = poc_so(idx_saz);
    poc_pfz = poc_so(idx_pfz);
    poc_az = poc_so(idx_az);
    eval(['poc_so_' all_abbs{i} '=poc_so(:);']);
    eval(['poc_stz_' all_abbs{i} '=poc_stz;']);
    eval(['poc_saz_' all_abbs{i} '=poc_saz;']);
    eval(['poc_pfz_' all_abbs{i} '=poc_pfz;']);
    eval(['poc_az_' all_abbs{i} '=poc_az;']);
end

rglist = {'so','stz','saz','pfz','az'};

valid_idx_so = ~isnan(poc_so_obs) & ~isnan(poc_so_cmc) & ~isnan(poc_so_can) & ~isnan(poc_so_gfd) & ...
    ~isnan(poc_so_ips) & ~isnan(poc_so_mhm) & ~isnan(poc_so_mhr) & ~isnan(poc_so_mlr) & ...
    ~isnan(poc_so_nlm) & ~isnan(poc_so_nmm);

valid_idx_stz = ~isnan(poc_stz_obs) & ~isnan(poc_stz_cmc) & ~isnan(poc_stz_can) & ~isnan(poc_stz_gfd) & ...
    ~isnan(poc_stz_ips) & ~isnan(poc_stz_mhm) & ~isnan(poc_stz_mhr) & ~isnan(poc_stz_mlr) & ...
    ~isnan(poc_stz_nlm) & ~isnan(poc_stz_nmm);

valid_idx_saz = ~isnan(poc_saz_obs) & ~isnan(poc_saz_cmc) & ~isnan(poc_saz_can) & ~isnan(poc_saz_gfd) & ...
    ~isnan(poc_saz_ips) & ~isnan(poc_saz_mhm) & ~isnan(poc_saz_mhr) & ~isnan(poc_saz_mlr) & ...
    ~isnan(poc_saz_nlm) & ~isnan(poc_saz_nmm);

valid_idx_pfz = ~isnan(poc_pfz_obs) & ~isnan(poc_pfz_cmc) & ~isnan(poc_pfz_can) & ~isnan(poc_pfz_gfd) & ...
    ~isnan(poc_pfz_ips) & ~isnan(poc_pfz_mhm) & ~isnan(poc_pfz_mhr) & ~isnan(poc_pfz_mlr) & ...
    ~isnan(poc_pfz_nlm) & ~isnan(poc_pfz_nmm);

valid_idx_az = ~isnan(poc_az_obs) & ~isnan(poc_az_cmc) & ~isnan(poc_az_can) & ~isnan(poc_az_gfd) & ...
    ~isnan(poc_az_ips) & ~isnan(poc_az_mhm) & ~isnan(poc_az_mhr) & ~isnan(poc_az_mlr) & ...
    ~isnan(poc_az_nlm) & ~isnan(poc_az_nmm);

for j = 1:length(rglist)
    for i = 1:length(all_abbs)
        eval(['poc_' rglist{j} '_' all_abbs{i} '_nonan=poc_' rglist{j} '_' all_abbs{i} '(valid_idx_' rglist{j} ');']);
    end
end

stats = [];
for j = 1:length(rglist)
    for i = 1:length(model_abbs)
        S = SStats(eval(['poc_' rglist{j} '_obs_nonan']),eval(['poc_' rglist{j} '_' model_abbs{i} '_nonan']));
        stats(j,:,i+1) = S;
    end
    stats(j,:,1) = SStats(eval(['poc_' rglist{j} '_obs_nonan']),eval(['poc_' rglist{j} '_obs_nonan']));

    stdev = stats(j,2,1); 

    stats(j,2,:) = stats(j,2,:)/stdev; 
    stats(j,3,:) = stats(j,3,:)/stdev;
end

stats(:,1,:) = [];

% Taylor diagram for surface POC (Fig.S5)
td_longname = {' ','ACCESS-ESM1-5','CanESM5','CESM2','CMCC-ESM2','CNRM-ESM2-1','GFDL-ESM4',...
    'IPSL-CM6A-LR','MIROC-ES2L','MPI-ESM-1-2-HAM','MPI-ESM1-2-HR',...
    'MPI-ESM1-2-LR','NorESM2-LM','NorESM2-MM','UKESM1-0-LL'};

figs5=figure(3);
set(gcf,'Position',[100 100 1500 1000])
subplot(2,3,1)
[hp, ht, axl] = taylor_diagram(squeeze(stats(1,1,:)),squeeze(stats(1,2,:)),squeeze(stats(1,3,:)), ...
    'markerLabel',td_longname,'markerLegend', 'on', ...
    'styleSTD', '-', 'colOBS','k', 'markerObs','o', ...
    'markerSize',6, 'tickRMS',0:1:2, ...
    'tickRMSangle', 115, 'showlabelsRMS', 'on', ...
    'titleRMS','on', 'titleOBS','Obs');
t1=title('Surface POC in SO','FontSize',14);
t1.Position(2) = t1.Position(2) + 1;

subplot(2,3,2)
[hp, ht, axl] = taylor_diagram(squeeze(stats(2,1,:)),squeeze(stats(2,2,:)),squeeze(stats(2,3,:)), ...
    'markerLabel',td_longname, 'markerLegend', 'on', ...
    'styleSTD', '-', 'colOBS','k', 'markerObs','o', ...
    'markerSize',6, 'tickRMS',0:1:3, ...
    'tickRMSangle', 115, 'showlabelsRMS', 'on', ...
    'titleRMS','on', 'titleOBS','Obs');
t2=title('Surface POC in STZ','FontSize',14);
t2.Position(2) = t2.Position(2) + 1.5;

subplot(2,3,3)
[hp, ht, axl] = taylor_diagram(squeeze(stats(3,1,:)),squeeze(stats(3,2,:)),squeeze(stats(3,3,:)), ...
    'markerLabel',td_longname, 'markerLegend', 'on', ...
    'styleSTD', '-', 'colOBS','k', 'markerObs','o', ...
    'markerSize',6, 'tickRMS',0:2:4, ...
    'tickRMSangle', 115, 'showlabelsRMS', 'on', ...
    'titleRMS','on', 'titleOBS','Obs');
t3=title('Surface POC in SAZ','FontSize',14);
t3.Position(2) = t3.Position(2) + 2;

subplot(2,3,4)
[hp, ht, axl] = taylor_diagram(squeeze(stats(4,1,:)),squeeze(stats(4,2,:)),squeeze(stats(4,3,:)), ...
    'markerLabel',td_longname, 'markerLegend', 'on', ...
    'styleSTD', '-', 'colOBS','k', 'markerObs','o', ...
    'markerSize',6, 'tickRMS',0:1:3, ...
    'tickRMSangle', 115, 'showlabelsRMS', 'on', ...
    'titleRMS','on', 'titleOBS','Obs');
t4=title('Surface POC in PFZ','FontSize',14);
t4.Position(2) = t4.Position(2) + 1.5;

subplot(2,3,5)
[hp, ht, axl] = taylor_diagram(squeeze(stats(5,1,:)),squeeze(stats(5,2,:)),squeeze(stats(5,3,:)), ...
    'markerLabel',td_longname, 'markerLegend', 'on', ...
    'styleSTD', '-', 'colOBS','k', 'markerObs','o', ...
    'markerSize',6, 'tickRMS',0:1:2, ...
    'tickRMSangle', 115, 'showlabelsRMS', 'on', ...
    'titleRMS','on', 'titleOBS','Obs');
t5=title('Surface POC in AZ','FontSize',14);
t5.Position(2) = t5.Position(2) + 1;

print(figs5,'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Figure\figs5','-djpeg','-r300');

% Calculate mean profiles for each zone
for i = 1:length(all_abbs)
    poc_djj_mean = eval(['poc_' all_abbs{i} '_djj_mean']);
    for k = 1:length(lev)
        poc_dep = poc_djj_mean(:,:,k);
        pocp_so(i,k) = nanmean(poc_dep(:));
        pocp_stz(i,k) = nanmean(poc_dep(idx_stz));
        pocp_saz(i,k) = nanmean(poc_dep(idx_saz));
        pocp_pfz(i,k) = nanmean(poc_dep(idx_pfz));
        pocp_az(i,k) = nanmean(poc_dep(idx_az));
    end
end

% Calculate threshold depth for poc
for j = 1:length(rglist)
    pocp_mat = eval(['pocp_' rglist{j}]);
    for i = 1:length(all_abbs)
        pocp_array = pocp_mat(i,:);
        idx = find(pocp_array < max(pocp_array)*0.5,1,'first');
        lev_thr = interp1([pocp_array(idx-1);pocp_array(idx)],[lev(idx-1);lev(idx)],max(pocp_array)*0.5);
        lev_thr_pocp(i,j) = lev_thr;
    end
end

% POC mean profile (Fig.S7)
splist2=[1;2;3;4;5;7;8;9;10;11;13;14;15;16;17];
figs7=figure(4);
set(gcf,'Position',[100 100 1200 800])
for s = 1:length(all_abbs)
    sp=subplot(3,6,splist2(s));
    plot(pocp_so(s,:),lev,'k-','LineWidth',2)
    hold on;
    plot(pocp_stz(s,:),lev,'m-','LineWidth',2)
    plot(pocp_saz(s,:),lev,'g-','LineWidth',2)
    plot(pocp_pfz(s,:),lev,'r-','LineWidth',2)
    plot(pocp_az(s,:),lev,'b-','LineWidth',2)
    xl=xlim;
    plot(xl,[lev_thr_pocp(s,1),lev_thr_pocp(s,1)],'k--','LineWidth',.5)
    plot(xl,[lev_thr_pocp(s,2),lev_thr_pocp(s,2)],'m--','LineWidth',.5)
    plot(xl,[lev_thr_pocp(s,3),lev_thr_pocp(s,3)],'g--','LineWidth',.5)
    plot(xl,[lev_thr_pocp(s,4),lev_thr_pocp(s,4)],'r--','LineWidth',.5)
    plot(xl,[lev_thr_pocp(s,5),lev_thr_pocp(s,5)],'b--','LineWidth',.5)
    hold off;
    title(all_longname{s},'FontSize',12)
    box on;
    grid on;
    % xlim([0 40])
    ylim([0 200])
    set(gca,'XAxisLocation','top','YDir','reverse')
end

legend({'SO','STZ','SAZ','PFZ','AZ'},'Position',[0.8, 0.10, 0.1, 0.2],'FontSize',12)

annotation('textbox', [0.13 0.45 0.05 0.1], 'String', 'Depth(m)', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
    'EdgeColor', 'none', 'FontSize', 12,'Rotation',90)

annotation('textbox', [0.45 0.05 0.1 0.05], 'String', 'POC(mg/m^{3})', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
    'EdgeColor', 'none', 'FontSize', 12)

print(figs7,'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Figure\figs7','-djpeg','-r300');

% Calculate mean profiles for chl
for i = 1:length(all_abbs)
    chl_djj_mean = eval(['chl_' all_abbs{i} '_djj_mean']);
    for k = 1:length(lev)
        chl_dep = chl_djj_mean(:,:,k);
        chlp_so(i,k) = nanmean(chl_dep(:));
        chlp_stz(i,k) = nanmean(chl_dep(idx_stz));
        chlp_saz(i,k) = nanmean(chl_dep(idx_saz));
        chlp_pfz(i,k) = nanmean(chl_dep(idx_pfz));
        chlp_az(i,k) = nanmean(chl_dep(idx_az));
    end
end

% Calculate threshold depth for chl
for j = 1:length(rglist)
    chlp_mat = eval(['chlp_' rglist{j}]);
    for i = 1:length(all_abbs)
        chlp_array = chlp_mat(i,:);
        idx = find(chlp_array < max(chlp_array)*0.1,1,'first');
        lev_thr = interp1([chlp_array(idx-1);chlp_array(idx)],[lev(idx-1);lev(idx)],max(chlp_array)*0.1);
        lev_thr_chlp(i,j) = lev_thr;
    end
end

% plot mean chlorophyll profile (Fig.S6)
figs6=figure(5);
set(gcf,'Position',[100 100 1200 800])
for s = 1:length(all_abbs)
    sp=subplot(3,6,splist2(s));
    plot(chlp_so(s,:),lev,'k-','LineWidth',2)
    hold on;
    plot(chlp_stz(s,:),lev,'m-','LineWidth',2)
    plot(chlp_saz(s,:),lev,'g-','LineWidth',2)
    plot(chlp_pfz(s,:),lev,'r-','LineWidth',2)
    plot(chlp_az(s,:),lev,'b-','LineWidth',2)
    xl=xlim;
    plot(xl,[lev_thr_chlp(s,1),lev_thr_chlp(s,1)],'k--','LineWidth',.5)
    plot(xl,[lev_thr_chlp(s,2),lev_thr_chlp(s,2)],'m--','LineWidth',.5)
    plot(xl,[lev_thr_chlp(s,3),lev_thr_chlp(s,3)],'g--','LineWidth',.5)
    plot(xl,[lev_thr_chlp(s,4),lev_thr_chlp(s,4)],'r--','LineWidth',.5)
    plot(xl,[lev_thr_chlp(s,5),lev_thr_chlp(s,5)],'b--','LineWidth',.5)
    hold off;
    title(all_longname{s},'FontSize',12)
    box on;
    grid on;
    % xlim([0 40])
    ylim([0 200])
    set(gca,'XAxisLocation','top','YDir','reverse')
end

legend({'SO','STZ','SAZ','PFZ','AZ'},'Position',[0.8, 0.1, 0.1, 0.2],'FontSize',12)

annotation('textbox', [0.13 0.45 0.05 0.1], 'String', 'Depth(m)', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
    'EdgeColor', 'none', 'FontSize', 12,'Rotation',90)

annotation('textbox', [0.45 0.05 0.1 0.05], 'String', 'Chlorophyll(mg/m^{3})', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
    'EdgeColor', 'none', 'FontSize', 12)

print(figs6,'D:\OneDrive - Australian National University\PhD\Project2\Model comparison\Figure\figs6','-djpeg','-r300');

% Calculate mean integrated phyc, zooc, and detoc in top 100m (Tab.S6-8)
for i = 1:length(all_abbs)
    chl_djj_mean = eval(['chl_' all_abbs{i} '_djj_mean']);
    chl_djj_mean2 = NaN(size(chl_djj_mean,1),size(chl_djj_mean,2),length(lev2));
    chl_djj_mean2(:,:,1) = chl_djj_mean(:,:,1);
    chl_djj_mean2(:,:,2:end) = chl_djj_mean;
    chl_djj_mean2(isnan(chl_djj_mean2)) = 0;
    chlz = NaN(size(chl_djj_mean,1),size(chl_djj_mean,2));
    for ji = 1:size(chl_djj_mean,1)
        for jj = 1:size(chl_djj_mean,2)
            chlz(ji,jj) = trapz(lev2(1:18),chl_djj_mean2(ji,jj,1:18));
        end
    end
    chlz(chlz==0) = NaN;
    mean_chl(i,1)=sum(chlz(~isnan(chlz)).*cosd(lat2r(~isnan(chlz))))...
        ./sum((cosd(lat2r(~isnan(chlz)))));
    mean_chl(i,2)=sum(chlz(~isnan(chlz)&idx_stz).*cosd(lat2r(~isnan(chlz)&idx_stz)))...
        ./sum(cosd(lat2r(~isnan(chlz)&idx_stz)));
    mean_chl(i,3)=sum(chlz(~isnan(chlz)&idx_saz).*cosd(lat2r(~isnan(chlz)&idx_saz)))...
        ./sum(cosd(lat2r(~isnan(chlz)&idx_saz)));
    mean_chl(i,4)=sum(chlz(~isnan(chlz)&idx_pfz).*cosd(lat2r(~isnan(chlz)&idx_pfz)))...
        ./sum(cosd(lat2r(~isnan(chlz)&idx_pfz)));
    mean_chl(i,5)=sum(chlz(~isnan(chlz)&idx_az).*cosd(lat2r(~isnan(chlz)&idx_az)))...
        ./sum(cosd(lat2r(~isnan(chlz)&idx_az)));
end

for i = 1:length(model_abbs)
    phyc_djj_mean = eval(['phyc_' model_abbs{i} '_djj_mean']);
    phyc_djj_mean2 = NaN(size(phyc_djj_mean,1),size(phyc_djj_mean,2),length(lev2));
    phyc_djj_mean2(:,:,1) = phyc_djj_mean(:,:,1);
    phyc_djj_mean2(:,:,2:end) = phyc_djj_mean;
    phyc_djj_mean2(isnan(phyc_djj_mean2)) = 0;
    phycz = NaN(size(phyc_djj_mean,1),size(phyc_djj_mean,2));
    for ji = 1:size(phyc_djj_mean,1)
        for jj = 1:size(phyc_djj_mean,2)
            phycz(ji,jj) = trapz(lev2(1:18),phyc_djj_mean2(ji,jj,1:18));
        end
    end
    phycz(phycz==0) = NaN;
    mean_phyc(i,1)=sum(phycz(~isnan(phycz)).*cosd(lat2r(~isnan(phycz))))...
        ./sum((cosd(lat2r(~isnan(phycz)))));
    mean_phyc(i,2)=sum(phycz(~isnan(phycz)&idx_stz).*cosd(lat2r(~isnan(phycz)&idx_stz)))...
        ./sum(cosd(lat2r(~isnan(phycz)&idx_stz)));
    mean_phyc(i,3)=sum(phycz(~isnan(phycz)&idx_saz).*cosd(lat2r(~isnan(phycz)&idx_saz)))...
        ./sum(cosd(lat2r(~isnan(phycz)&idx_saz)));
    mean_phyc(i,4)=sum(phycz(~isnan(phycz)&idx_pfz).*cosd(lat2r(~isnan(phycz)&idx_pfz)))...
        ./sum(cosd(lat2r(~isnan(phycz)&idx_pfz)));
    mean_phyc(i,5)=sum(phycz(~isnan(phycz)&idx_az).*cosd(lat2r(~isnan(phycz)&idx_az)))...
        ./sum(cosd(lat2r(~isnan(phycz)&idx_az)));
end

for i = 1:length(model_abbs)
    zooc_djj_mean = eval(['zooc_' model_abbs{i} '_djj_mean']);
    zooc_djj_mean2 = NaN(size(zooc_djj_mean,1),size(zooc_djj_mean,2),length(lev2));
    zooc_djj_mean2(:,:,1) = zooc_djj_mean(:,:,1);
    zooc_djj_mean2(:,:,2:end) = zooc_djj_mean;
    zooc_djj_mean2(isnan(zooc_djj_mean2)) = 0;
    zoocz = NaN(size(zooc_djj_mean,1),size(zooc_djj_mean,2));
    for ji = 1:size(zooc_djj_mean,1)
        for jj = 1:size(zooc_djj_mean,2)
            zoocz(ji,jj) = trapz(lev2(1:18),zooc_djj_mean2(ji,jj,1:18));
        end
    end
    zoocz(zoocz==0) = NaN;
    mean_zooc(i,1)=sum(zoocz(~isnan(zoocz)).*cosd(lat2r(~isnan(zoocz))))...
        ./sum((cosd(lat2r(~isnan(zoocz)))));
    mean_zooc(i,2)=sum(zoocz(~isnan(zoocz)&idx_stz).*cosd(lat2r(~isnan(zoocz)&idx_stz)))...
        ./sum(cosd(lat2r(~isnan(zoocz)&idx_stz)));
    mean_zooc(i,3)=sum(zoocz(~isnan(zoocz)&idx_saz).*cosd(lat2r(~isnan(zoocz)&idx_saz)))...
        ./sum(cosd(lat2r(~isnan(zoocz)&idx_saz)));
    mean_zooc(i,4)=sum(zoocz(~isnan(zoocz)&idx_pfz).*cosd(lat2r(~isnan(zoocz)&idx_pfz)))...
        ./sum(cosd(lat2r(~isnan(zoocz)&idx_pfz)));
    mean_zooc(i,5)=sum(zoocz(~isnan(zoocz)&idx_az).*cosd(lat2r(~isnan(zoocz)&idx_az)))...
        ./sum(cosd(lat2r(~isnan(zoocz)&idx_az)));
end

for i = 1:length(model_abbs)
    detoc_djj_mean = eval(['detoc_' model_abbs{i} '_djj_mean']);
    detoc_djj_mean2 = NaN(size(detoc_djj_mean,1),size(detoc_djj_mean,2),length(lev2));
    detoc_djj_mean2(:,:,1) = detoc_djj_mean(:,:,1);
    detoc_djj_mean2(:,:,2:end) = detoc_djj_mean;
    detoc_djj_mean2(isnan(detoc_djj_mean2)) = 0;
    detocz = NaN(size(detoc_djj_mean,1),size(detoc_djj_mean,2));
    for ji = 1:size(detoc_djj_mean,1)
        for jj = 1:size(detoc_djj_mean,2)
            detocz(ji,jj) = trapz(lev2(1:18),detoc_djj_mean2(ji,jj,1:18));
        end
    end
    detocz(detocz==0) = NaN;
    mean_detoc(i,1)=sum(detocz(~isnan(detocz)).*cosd(lat2r(~isnan(detocz))))...
        ./sum((cosd(lat2r(~isnan(detocz)))));
    mean_detoc(i,2)=sum(detocz(~isnan(detocz)&idx_stz).*cosd(lat2r(~isnan(detocz)&idx_stz)))...
        ./sum(cosd(lat2r(~isnan(detocz)&idx_stz)));
    mean_detoc(i,3)=sum(detocz(~isnan(detocz)&idx_saz).*cosd(lat2r(~isnan(detocz)&idx_saz)))...
        ./sum(cosd(lat2r(~isnan(detocz)&idx_saz)));
    mean_detoc(i,4)=sum(detocz(~isnan(detocz)&idx_pfz).*cosd(lat2r(~isnan(detocz)&idx_pfz)))...
        ./sum(cosd(lat2r(~isnan(detocz)&idx_pfz)));
    mean_detoc(i,5)=sum(detocz(~isnan(detocz)&idx_az).*cosd(lat2r(~isnan(detocz)&idx_az)))...
        ./sum(cosd(lat2r(~isnan(detocz)&idx_az)));
end

for i = 1:length(bac_abbs)
    bacc_djj_mean = eval(['bacc_' bac_abbs{i} '_djj_mean']);
    bacc_djj_mean2 = NaN(size(bacc_djj_mean,1),size(bacc_djj_mean,2),length(lev2));
    bacc_djj_mean2(:,:,1) = bacc_djj_mean(:,:,1);
    bacc_djj_mean2(:,:,2:end) = bacc_djj_mean;
    bacc_djj_mean2(isnan(bacc_djj_mean2)) = 0;
    baccz = NaN(size(bacc_djj_mean,1),size(bacc_djj_mean,2));
    for ji = 1:size(bacc_djj_mean,1)
        for jj = 1:size(bacc_djj_mean,2)
            baccz(ji,jj) = trapz(lev2(1:18),bacc_djj_mean2(ji,jj,1:18));
        end
    end
    baccz(baccz==0) = NaN;
    mean_bacc(i,1)=sum(baccz(~isnan(baccz)).*cosd(lat2r(~isnan(baccz))))...
        ./sum((cosd(lat2r(~isnan(baccz)))));
    mean_bacc(i,2)=sum(baccz(~isnan(baccz)&idx_stz).*cosd(lat2r(~isnan(baccz)&idx_stz)))...
        ./sum(cosd(lat2r(~isnan(baccz)&idx_stz)));
    mean_bacc(i,3)=sum(baccz(~isnan(baccz)&idx_saz).*cosd(lat2r(~isnan(baccz)&idx_saz)))...
        ./sum(cosd(lat2r(~isnan(baccz)&idx_saz)));
    mean_bacc(i,4)=sum(baccz(~isnan(baccz)&idx_pfz).*cosd(lat2r(~isnan(baccz)&idx_pfz)))...
        ./sum(cosd(lat2r(~isnan(baccz)&idx_pfz)));
    mean_bacc(i,5)=sum(baccz(~isnan(baccz)&idx_az).*cosd(lat2r(~isnan(baccz)&idx_az)))...
        ./sum(cosd(lat2r(~isnan(baccz)&idx_az)));
end

for i = 1:length(all_abbs)
    poc_djj_mean = eval(['poc_' all_abbs{i} '_djj_mean']);
    poc_djj_mean2 = NaN(size(poc_djj_mean,1),size(poc_djj_mean,2),length(lev2));
    poc_djj_mean2(:,:,1) = poc_djj_mean(:,:,1);
    poc_djj_mean2(:,:,2:end) = poc_djj_mean;
    poc_djj_mean2(isnan(poc_djj_mean2)) = 0;
    pocz = NaN(size(poc_djj_mean,1),size(poc_djj_mean,2));
    for ji = 1:size(poc_djj_mean,1)
        for jj = 1:size(poc_djj_mean,2)
            pocz(ji,jj) = trapz(lev2(1:18),poc_djj_mean2(ji,jj,1:18));
        end
    end
    pocz(pocz==0) = NaN;
    mean_poc(i,1)=sum(pocz(~isnan(pocz)).*cosd(lat2r(~isnan(pocz))))...
        ./sum((cosd(lat2r(~isnan(pocz)))));
    mean_poc(i,2)=sum(pocz(~isnan(pocz)&idx_stz).*cosd(lat2r(~isnan(pocz)&idx_stz)))...
        ./sum(cosd(lat2r(~isnan(pocz)&idx_stz)));
    mean_poc(i,3)=sum(pocz(~isnan(pocz)&idx_saz).*cosd(lat2r(~isnan(pocz)&idx_saz)))...
        ./sum(cosd(lat2r(~isnan(pocz)&idx_saz)));
    mean_poc(i,4)=sum(pocz(~isnan(pocz)&idx_pfz).*cosd(lat2r(~isnan(pocz)&idx_pfz)))...
        ./sum(cosd(lat2r(~isnan(pocz)&idx_pfz)));
    mean_poc(i,5)=sum(pocz(~isnan(pocz)&idx_az).*cosd(lat2r(~isnan(pocz)&idx_az)))...
        ./sum(cosd(lat2r(~isnan(pocz)&idx_az)));
end

