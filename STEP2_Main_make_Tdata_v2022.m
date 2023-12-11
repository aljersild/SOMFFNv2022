% To create the multidimensional dataset from specific datasets, one for
% each parameter. Datasets for each parameter should be rescaled by this
% point. Training dataset structure must be the same as labelling dataset.
% This code follows the work of Maciek Telszewski and loads SST, MLD, CHL, 
% BATHY, etc. data from e.g. reanalysis products to train the NN 

%Peter Landschutzer 17.01.2012
%University of East Anglia, Norwich
%Edited by Annika Jersild July 2022

Input_Training_and_Labelling_new_GUI_v2022
% -------------------------------------------------------------------------
% 1. define if global or Atlantic data are to prepare using a GUI
% -------------------------------------------------------------------------

prompt={'ocean'};    
defans={'please enter "Atlantic" or "Global"'};
fields = {'ocean'};  
info = inputdlg(prompt, 'please enter if global or Atlantic', 1, defans);
if ~isempty(info)              
   info = cell2struct(info,fields);
   ocean    = info.ocean;
end      
clear prompt defans fields info

% -------------------------------------------------------------------------
% 2. load the data which are in (264,180,360) format
% -------------------------------------------------------------------------

% V2 used: ECCO, Raynolds, Globalviev, Taka, Globecolor

load('training_data/atm_co2/atm_co2_grid_v2022.mat');
atm_co2=data_all;
clear data_all lat lon
load('training_data/mld/mld_clim_v2022.mat');
clear lonmld latmld
load('training_data/salinity/sss_v2022.mat');
clear lonsss latsss
load('training_data/taka/Taka_pCO2_eth_v2022.mat');
load('training_data/sst/sst_v2022.mat');
load('training_data/chl/chl_v2022.mat');
load('output/BIOMEoutput_SOCAT/networks/SOM_biome_4x4.mat');


for year2go=min(year_output):max(year_output)

    year2go
    
% -------------------------------------------------------------------------
% 3. The first parameter is SST
% -------------------------------------------------------------------------

%========================= work model sst =================================

y_ind=find(years==year2go);
sst1=sst(y_ind,:,:);
time1=time(y_ind,:,:);
%months1=str2double(datestr(time1,'mm'));
latsst1=latsst(y_ind,:,:);
lonsst1=lonsst(y_ind,:,:);



%%%%%%%%%%%%%%%%%%%%%
%cropping
%%%%%%%%%%%%%%%%%%%%%
%load('lsmask.mat');
%dx=d1;
%dx(:,1:end/2)=d1(:,end/2+1:end);
%dx(:,end/2+1:end)=d1(:,1:end/2);
%lon=lonsst1;
%lat=latsst1;
%sst1(:,dx<invareamin)=NaN;
%sst1(:,dx>invareamax)=NaN;

%[sst1]=crop_field(sst1,invareamin,invareamax);

%clear d1

[sst1, lon1, lat1]=latlon_crop(sst1,lon,lat,loncropmin, loncropmax,latcropmin,latcropmax);
[time1, lon1, lat1]=latlon_crop(time1,lon,lat,loncropmin, loncropmax,latcropmin,latcropmax);
[latsst1, lon1, lat1]=latlon_crop(latsst1,lon,lat,loncropmin, loncropmax,latcropmin,latcropmax);
[lonsst1, lon1, lat1]=latlon_crop(lonsst1,lon,lat,loncropmin, loncropmax,latcropmin,latcropmax);

%=============== sst is 3 dim array - reshape to 1D =======================

sst_new = reshape(sst1, prod(size(sst1)), 1);
time_new = reshape(time1, prod(size(time1)), 1);
lat_new = reshape(latsst1, prod(size(latsst1)), 1);
lon_new = reshape(lonsst1, prod(size(lonsst1)), 1);

%============== create an array with the same length for year =============

for i=1:1:length(sst_new);
    year2go_new(i)=year2go;
end
year2go_new=year2go_new';

%=================== put data in the same array ===========================
% 1=year; 2=matlab time 3=lat; 4=lon; 5 =sst

sst_NA_ALL=[year2go_new time_new lat_new lon_new sst_new];

%========= Insert 1 column in front of the existing data for months: ======

sst_NA_ALL = sst_NA_ALL(:,[1 1 2 3 4 5]); 

% column 1 for year

% column 2 for month:

months=zeros(12,180,360);
months(:,:,:)=NaN;
for uu=1:12
    months(uu,1:end,1:end)=uu;
end

[months1, lon1, lat1]=latlon_crop(months,lon,lat,loncropmin, loncropmax,latcropmin,latcropmax);
months_new = reshape(months1, prod(size(months1)), 1);
sst_NA_ALL(:,2)=months_new;

clear time1 i months;

%==== Insert columns at the end of the existing data for parameter: =======
% 7=mld 8=chl 9=sss 10=atm_co2 11=data_taka 12=sst anom
% 13=mld anom 14 =chl anom 15 =sss anom 16=atm_co2 anom 17=NaN 18=biomes

sst_NA_ALL = sst_NA_ALL(:,[1 2 3 4 5 6 6 6 6 6 6 6 6 6 6 6 6 6]); 

sst_NA_ALL(:,7:18)=NaN;
% sst_NA_ALL(:,8)=NaN;
% sst_NA_ALL(:,9)=NaN;
% sst_NA_ALL(:,10)=NaN;
% sst_NA_ALL(:,11)=NaN;
% sst_NA_ALL(:,12)=NaN;
% sst_NA_ALL(:,13)=NaN;
% sst_NA_ALL(:,14)=NaN;
% sst_NA_ALL(:,15)=NaN;
% sst_NA_ALL(:,16)=NaN;
% sst_NA_ALL(:,17)=NaN;
% sst_NA_ALL(:,18)=NaN;
% sst_NA_ALL(:,19)=NaN;
% sst_NA_ALL(:,20)=NaN;
% sst_NA_ALL(:,21)=NaN;
% sst_NA_ALL(:,22)=NaN;
% sst_NA_ALL(:,23)=NaN;
% sst_NA_ALL(:,24)=NaN;

%--------------------------------------------------------------------------
% 4. Coinciding training data
%--------------------------------------------------------------------------

% bathy1=bathy(y_ind,:,:);
mld1=mld(y_ind,:,:);
atm_co21=atm_co2(y_ind,:,:);
% spec_rate1=spec_rate(y_ind,:,:);
chl1=chl(y_ind,:,:);
sss1=sss(y_ind,:,:);
data_taka1=data_taka(y_ind,:,:);
% wind1=wind(y_ind,:,:);
%ssh1=ssh(y_ind,:,:);
%curr_speed1=curr_speed(y_ind,:,:);
%mld_anom1=mld_anom(y_ind,:,:);
sst_anom1=sst_anom(y_ind,:,:);
chl_anom1=chl_anom(y_ind,:,:);
sss_anom1=sss_anom(y_ind,:,:);
atm_anom1=data_anom(y_ind,:,:);
biomes1=biomes(y_ind,:,:);
% 
% % 
% % bathy1(:,dx<invareamin)=NaN;
% % bathy1(:,dx>invareamax)=NaN;
% mld1(:,dx<invareamin)=NaN;
% mld1(:,dx>invareamax)=NaN;
% atm_co21(:,dx<invareamin)=NaN;
% atm_co21(:,dx>invareamax)=NaN;
% % spec_rate1(:,dx<invareamin)=NaN;
% % spec_rate1(:,dx>invareamax)=NaN;
% chl1(:,dx<invareamin)=NaN;
% chl1(:,dx>invareamax)=NaN;
% sss1(:,dx<invareamin)=NaN;
% sss1(:,dx>invareamax)=NaN;
% data_taka1(:,dx<invareamin)=NaN;
% data_taka1(:,dx>invareamax)=NaN;
% % wind1(:,dx<invareamin)=NaN;
% % wind1(:,dx>invareamax)=NaN;
% % ssh1(:,dx<invareamin)=NaN;
% % ssh1(:,dx>invareamax)=NaN;
% % curr_speed1(:,dx<invareamin)=NaN;
% % curr_speed1(:,dx>invareamax)=NaN;
% %mld_anom1(:,dx<invareamin)=NaN;
% %mld_anom1(:,dx>invareamax)=NaN;
% chl_anom1(:,dx<invareamin)=NaN;
% chl_anom1(:,dx>invareamax)=NaN;
% sst_anom1(:,dx<invareamin)=NaN;
% sst_anom1(:,dx>invareamax)=NaN;
% sss_anom1(:,dx<invareamin)=NaN;
% sss_anom1(:,dx>invareamax)=NaN;
% atm_anom1(:,dx<invareamin)=NaN;
% atm_anom1(:,dx>invareamax)=NaN;
% biomes1(:,dx<invareamin)=NaN;
% biomes1(:,dx>invareamax)=NaN;

[mld1, lonmld1, latmld1]=latlon_crop(mld1,lon,lat,loncropmin, loncropmax,latcropmin,latcropmax);
[atm_co21, lonatm_co21, latatm_co21]=latlon_crop(atm_co21,lon,lat,loncropmin, loncropmax,latcropmin,latcropmax);
% [bathy1, lonbathy1, latbathy1]=latlon_crop(bathy1,lon,lat,loncropmin, loncropmax,latcropmin,latcropmax);
[data_taka1, londata_taka1, latdata_taka1]=latlon_crop(data_taka1,lon,lat,loncropmin, loncropmax,latcropmin,latcropmax);
[sss1, lonsss1, latsss1]=latlon_crop(sss1,lon,lat,loncropmin, loncropmax,latcropmin,latcropmax);
% [spec_rate1, lonspec_rate1, latspec_rate1]=latlon_crop(spec_rate1,lon,lat,loncropmin, loncropmax,latcropmin,latcropmax);
[chl1, lonchl1, latchl1]=latlon_crop(chl1,lon,lat,loncropmin, loncropmax,latcropmin,latcropmax);
%[ssh1, lonssh1, latssh1]=latlon_crop(ssh1,lon,lat,loncropmin, loncropmax,latcropmin,latcropmax);
%[curr_speed1, loncurr_speed1, latcurr_speed1]=latlon_crop(curr_speed1,lon,lat,loncropmin, loncropmax,latcropmin,latcropmax);
% [wind1, lonwind1, latwind1]=latlon_crop(wind1,lon,lat,loncropmin, loncropmax,latcropmin,latcropmax);
%[mld_anom1, lonmld1, latmld1]=latlon_crop(mld_anom1,lon,lat,loncropmin, loncropmax,latcropmin,latcropmax);
[sst_anom1, lonsst11, latsst11]=latlon_crop(sst_anom1,lon,lat,loncropmin, loncropmax,latcropmin,latcropmax);
[sss_anom1, lonsss1, latsss1]=latlon_crop(sss_anom1,lon,lat,loncropmin, loncropmax,latcropmin,latcropmax);
[chl_anom1, lonchl1, latchl1]=latlon_crop(chl_anom1,lon,lat,loncropmin, loncropmax,latcropmin,latcropmax);
[atm_anom1, lonatm_anom1, latatm_anom1]=latlon_crop(atm_anom1,lon,lat,loncropmin, loncropmax,latcropmin,latcropmax);
[biomes1, lonbio_anom1, latbio_anom1]=latlon_crop(biomes1,lon,lat,loncropmin, loncropmax,latcropmin,latcropmax);


% bathy_new = reshape(bathy1, prod(size(bathy1)), 1);
mld_new = reshape(mld1, prod(size(mld1)), 1);
atm_co2_new = reshape(atm_co21, prod(size(atm_co21)), 1);
% spec_rate_new = reshape(spec_rate1, prod(size(spec_rate1)), 1);
chl_new = reshape(chl1, prod(size(chl1)), 1);
sss_new = reshape(sss1, prod(size(sss1)), 1);
data_taka_new = reshape(data_taka1, prod(size(data_taka1)), 1);
%ssh_new = reshape(ssh1, prod(size(ssh1)), 1);
%curr_speed_new = reshape(curr_speed1, prod(size(curr_speed1)), 1);
% wind_new = reshape(wind1, prod(size(wind1)), 1);
%mld_anom_new = reshape(mld_anom1, prod(size(mld_anom1)), 1);
chl_anom_new = reshape(chl_anom1, prod(size(chl_anom1)), 1);
sst_anom_new = reshape(sst_anom1, prod(size(sst_anom1)), 1);
sss_anom_new = reshape(sss_anom1, prod(size(sss_anom1)), 1);
atm_anom_new = reshape(atm_anom1, prod(size(atm_anom1)), 1);
biomes_new = reshape(biomes1, prod(size(biomes1)), 1);

% sst_NA_ALL(:,7)=bathy_new;
sst_NA_ALL(:,7)=mld_new;
sst_NA_ALL(:,8)=chl_new;
sst_NA_ALL(:,9)=sss_new;
sst_NA_ALL(:,10)=atm_co2_new;
% sst_NA_ALL(:,12)=spec_rate_new;
sst_NA_ALL(:,11)=data_taka_new;
% sst_NA_ALL(:,15)=wind_new;
%sst_NA_ALL(:,16)=ssh_new;
%sst_NA_ALL(:,17)=curr_speed_new;
sst_NA_ALL(:,12)=sst_anom_new;
%sst_NA_ALL(:,13)=mld_anom_new;
%sst_NA_ALL(:,13)=NaN;
sst_NA_ALL(:,14)=chl_anom_new;
sst_NA_ALL(:,15)=sss_anom_new;
sst_NA_ALL(:,16)=atm_anom_new;
sst_NA_ALL(:,18)=biomes_new;

clear a b c d i j k lat_new lon_new time_new sst_new year2go_new
clear bathy_new mld_new chl_new data_taka_new sss_new spec_rate_new atm_co2_new ssh_new curr_speed_new wind_new

Tdata_500m=sst_NA_ALL;
clear sst_NA_ALL;
% -------------------------------------------------------------------------
% 5. The next parameter is the sinosoidal representation of time
% -------------------------------------------------------------------------

% for i=1:length(Tdata_500m)
%     Tdata_500m(i,13)=sin(((-3/4)*pi)+((Tdata_500m(i,2)-0.5)* 0.1667)*pi);
% end

%======== order the data (as L), first in time, then in lon and lat =======

Tdata_500m = sortrows(Tdata_500m,[3 4 5]);

% 7=mld 8=chl 9=sss 10=atm_co2 11=data_taka 12=sst anom
% 13=mld anom 14 =chl anom 15 =sss anom 16=atm_co2 anom 17=NaN
Tdata_500m=Tdata_500m(:,[1 2 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18]);
%==================== still from Maciek ===================================

% Remaining land pixels come from the SST data where 99 denotes land
% according to JMA. There are only 5 such disputable pixels for the entire
% study area so I simply remove them:
% ind2=find(Tdata_500m(:,6)>90);
% Tdata_500m(ind2,:)=[];
% 
% clear ind2;

% -------------------------------------------------------------------------
% 6. save output
% -------------------------------------------------------------------------
% 1=yr 2=mon 3=time 4=lat 5=lon 6=sst 7=bathy 8=mld 9=chl 10=sss 11=atm_co2
% 12=spec 13=sinosoidal time 14=data_taka 15=wind 16=ssh 17=curr 20 sst
% anom 21 mld anom 22=chl anom 23=sss anom 24=atmco2 anom

Tdata=Tdata_500m;

%eval(['Tdata_' num2str(year2go) '_500m = Tdata_500m;']);
eval(['Tdata_' num2str(year2go) ' = Tdata;']);

if(strcmp(ocean,'Global')==1)
saveloc=['output/Tdata/Tdata_' num2str(year2go) '_v2022.mat'];
elseif(strcmp(ocean,'Atlantic')==1)
%saveloc=['output/Tdata/Tdata_' num2str(year2go) '_anom_v3_ML_biomes_30yr.mat'];
end

savedat=['Tdata_' num2str(year2go)];

save(saveloc, savedat);

clear saveloc savedat Tdata_500m

end
% -------------------------------------------------------------------------
% 7. graphical check if ok
% -------------------------------------------------------------------------

xx=size(lat);
yy=size(lon);

s = ['Tdata_' num2str(year2go)];
plot_array=eval(s);
clear s
indxx=find(plot_array(:,2)==month2plot & plot_array(:,1)==year2go);

LON=plot_array(indxx,4);
LAT=plot_array(indxx,3);

[field]=vec_to_array1(plot_array(indxx,5:17),lon,lat,LON,LAT);

% for ii=1:1:yy(1)
%    latt(:,ii) = lat(:);   
% end
% for ii=1:1:xx(1)
%    lonn(ii,:) = lon(:);  
% end

map_plot_2016(lon, lat, squeeze(field(1,:,:)), 'equidistant', lonmin, lonmax, latmin, latmax, 1, 0, 30,2, 'SST',1,2)
% map_plot(lon, lat, squeeze(field(2,:,:)), 'equidistant', lonmin, lonmax, latmin, latmax, 2, -5000, 0, 'BATHY',1)
map_plot_2016(lon, lat, squeeze(field(2,:,:)), 'equidistant', lonmin, lonmax, latmin, latmax, 3, 0, 8,1, 'MLD',1,2)
map_plot_2016(lon, lat, squeeze(field(3,:,:)), 'equidistant', lonmin, lonmax, latmin, latmax, 4, -4, 4,1, 'CHL',1,2)
map_plot_2016(lon, lat, squeeze(field(4,:,:)), 'equidistant', lonmin, lonmax, latmin, latmax, 5, 32, 40,1, 'SSS',1,2)
map_plot_2016(lon, lat, squeeze(field(5,:,:)), 'equidistant', lonmin, lonmax, latmin, latmax, 6, 350, 420,2, 'ATM_CO2',1,2)
%map_plot(lon, lat, squeeze(field(7,:,:)), 'equidistant', lonmin, lonmax, latmin, latmax, 7, 0, 16, 'biomes',1)
map_plot_2016(lon, lat, squeeze(field(6,:,:)), 'equidistant', lonmin, lonmax, latmin, latmax, 8, 300, 450,20, 'data_taka',1,2)
%map_plot(lon, lat, squeeze(field(10,:,:)), 'equidistant', lonmin, lonmax, latmin, latmax, 9, 0, 100, 'wind',1)
%map_plot(lon, lat, squeeze(field(11,:,:)), 'equidistant', lonmin, lonmax, latmin, latmax, 10, -2, 2, 'ssh',1)
%map_plot(lon, lat, squeeze(field(12,:,:)), 'equidistant', lonmin, lonmax, latmin, latmax, 11, 0, 1, 'curr speed',1)
map_plot_2016(lon, lat, squeeze(field(7,:,:)), 'equidistant', lonmin, lonmax, latmin, latmax, 12, -2, 2,0.5, 'SST anom',1,2)
map_plot_new(lon, lat, squeeze(field(8,:,:)), 'equidistant', lonmin, lonmax, latmin, latmax, 13, -2, 2, 'mld anom',1,2)
map_plot_2016(lon, lat, squeeze(field(9,:,:)), 'equidistant', lonmin, lonmax, latmin, latmax, 14, -2, 2,0.5, 'chl anom',1,2)
map_plot_2016(lon, lat, squeeze(field(10,:,:)), 'equidistant', lonmin, lonmax, latmin, latmax, 15, -2, 2,0.5, 'sss anom',1,2)
map_plot_2016(lon, lat, squeeze(field(11,:,:)), 'equidistant', lonmin, lonmax, latmin, latmax, 16, 0, 44,0.2, 'atm anom',1,2)
map_plot_2016(lon, lat, squeeze(field(13,:,:)), 'equidistant', lonmin, lonmax, latmin, latmax, 17, 0, 50,2, 'biomes',1,2)