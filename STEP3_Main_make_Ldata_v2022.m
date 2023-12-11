% This is the main pre-NN in-situ data handling file.
% This code follows the work of Maciek Telszewski and creates for every
% year a summary file of the observations. Every column consists of one
% training parameter for the NN and additionally SOCAT fCO2 parameter from
% the SOCAT voyages. LON LAT YR and MONTH for the data are determined by
% input file

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
% 2. load observation files
% -------------------------------------------------------------------------

for year2go=min(year_output):max(year_output)

    year2go
    
if(strcmp(ocean,'Global')==1)    

path =(['observations/SOCATv2022/' num2str(year2go) '/data_' num2str(year2go) '_v2022.mat']);
%path =(['observations/SOCATv6/' num2str(year2go) '/data_' num2str(year2go) '_SOCATv6_v2018_SOCCOM_offset.mat']);

elseif(strcmp(ocean,'Atlantic')==1)
path =(['observations/SOCATv2.0/' num2str(year2go) '/data_' num2str(year2go) '_SOCAT.mat']);
%path =(['observations/SOCAT4GCB/' num2str(year2go) '/data_' num2str(year2go) '_SOCAT_global_pco2V2_GBC.mat']);
end
load(path)

% creating a temp file just in case sth goes wrong.
sk = ['data_' num2str(year2go) '_SOCAT'];
%sk = ['data_' num2str(year2go) '_SOCAT'];
temp=eval(sk);
clear sk path

% for Ldata_com
% if year2go==2011
%     load(path2)
%     sl=['data_' num2str(year2go) '_NIESUEA'];
%     tteemmpp=eval(sl);
%     temp=[temp;tteemmpp];
%     clear sl path1
% end

%======================= remove fco2 nans =================================

% column 5 is fCO2
ind = find(isnan(temp(:,5)));
temp(ind,:)=[]; clear ind;

%================= create new array with all parameters ===================

% columns: 1.year/2.month/3.lat/4.lon/5.fco2

Ldata_binned1=temp(:,[1 2 3 4 5]);
Ldata_binned1=double(Ldata_binned1);
%Ldata_binned1(:,3)=datenum(Ldata_binned1(:,1),Ldata_binned1(:,2),1);
Ldata_binned=unique(Ldata_binned1, 'rows');
clear temp Ldata_binned1;
% column 6 in Ldata_All is fco2
ind=find(Ldata_binned(:,5)<50);
Ldata_binned(ind,:)=[];

Ldata_binned(:,6)=NaN;
Ldata_coincided = Ldata_binned(:,[1 2 3 4 6 6 6 6 6 6 6 6 6 6 6 5 6]);

%======================= save in-between output ===========================

%save output/Ldata/Ldata_binned.mat Ldata_binned
%save output/Ldata/Ldata_All.mat Ldata_All

%clear Ldata_binned Ldata_All

% -------------------------------------------------------------------------
% 3. Coinciding other training parameter
% -------------------------------------------------------------------------

% for data incl. coast
if(strcmp(ocean,'Global')==1)
path =(['output/Tdata/Tdata_' num2str(year2go) '_v2022.mat']);
elseif(strcmp(ocean,'Atlantic')==1)
path =(['output/Tdata/Tdata_' num2str(year2go) '_anom_v4_ML_biomes_30yr.mat']);
end
load(path);
clear path
 
s = ['Tdata_' num2str(year2go)];
Tdata_500m=eval(s);

for i= 1:length(Ldata_coincided)

    [lat2]=coincide_LdataV2(Ldata_coincided,Tdata_500m,i);

    Ldata_coincided(i,5)=lat2(:,5);
    Ldata_coincided(i,6)=lat2(:,6);
    Ldata_coincided(i,7)=lat2(:,7);
    Ldata_coincided(i,8)=lat2(:,8);
    Ldata_coincided(i,9)=lat2(:,9);
    Ldata_coincided(i,10)=lat2(:,10);
    Ldata_coincided(i,11)=lat2(:,11);
    Ldata_coincided(i,12)=lat2(:,12);
    Ldata_coincided(i,13)=lat2(:,13);
    Ldata_coincided(i,14)=lat2(:,14);
    Ldata_coincided(i,15)=lat2(:,15);
    Ldata_coincided(i,17)=lat2(:,17);
end
clear lat1 lon1 mon1 mon2find mon2 lon2find lon2 lat2find lat2 i


clear Tdata_500m s

%======================= remove sst nans ==================================

ind2=find(isnan(Ldata_coincided(:,5)));
Ldata_coincided(ind2,:)=[];
 
clear ind2;

% -------------------------------------------------------------------------
% 4. save output
% -------------------------------------------------------------------------

eval(['Ldata_' num2str(year2go) '_coincided_v2 =  Ldata_coincided;']);

if(strcmp(ocean,'Global')==1)
%saveloc=['output/Ldata/Ldata_' num2str(year2go) '_coincided_v2017.mat'];
%saveloc=['output/Ldata/Ldata_' num2str(year2go) '_coincided_SOCATv6_v2018_SOCCOM_offset.mat'];
saveloc=['output/Ldata/Ldata_' num2str(year2go) '_v2022.mat'];
elseif(strcmp(ocean,'Atlantic')==1)
saveloc=['output/Ldata/Ldata_' num2str(year2go) '_coincided_v2_anom_v3_ML_biomes_30yr'];
end
savedat=['Ldata_' num2str(year2go) '_coincided_v2'];

save(saveloc, savedat);

clear saveloc savedat

end
% -------------------------------------------------------------------------
% 5. graphical check if ok
% -------------------------------------------------------------------------

lat =latcropmin+0.5:1:latcropmax;
lon =loncropmin+0.5:1:loncropmax;

xx=size(lat);
yy=size(lon);

month2plot=4;

indxx=find(Ldata_coincided(:,2)==month2plot & Ldata_coincided(:,1)==year2go);

LON=Ldata_coincided(indxx,4);
LAT=Ldata_coincided(indxx,3);

[field]=vec_to_array1(Ldata_coincided(indxx,5:17),lon,lat,LON,LAT);

% latt=zeros(xx(2),yy(2));
% lonn=zeros(xx(2),yy(2));
% for ii=1:1:yy(2)
%    latt(:,ii) = lat(:);   
% end
% for ii=1:1:xx(2)
%    lonn(ii,:) = lon(:);  
% end

map_plot_new(lon, lat, squeeze(field(1,:,:)), 'equidistant', lonmin, lonmax, latmin, latmax, 1, 0, 30, 'SST',1,2)
%map_plot(lon, lat, squeeze(field(2,:,:)), 'equidistant', lonmin, lonmax, latmin, latmax, 2, -5000, 0, 'BATHY',1)
map_plot_new(lon, lat, squeeze(field(2,:,:)), 'equidistant', lonmin, lonmax, latmin, latmax, 3, 0, 8, 'MLD',1,2)
map_plot_new(lon, lat, squeeze(field(3,:,:)), 'equidistant', lonmin, lonmax, latmin, latmax, 4, -4, 4, 'CHL',1,2)
map_plot_new(lon, lat, squeeze(field(4,:,:)), 'equidistant', lonmin, lonmax, latmin, latmax, 5, 30, 40, 'SSS',1,2)
map_plot_new(lon, lat, squeeze(field(5,:,:)), 'equidistant', lonmin, lonmax, latmin, latmax, 6, 350, 420, 'ATM_CO2',1,2)
%map_plot(lon, lat, squeeze(field(7,:,:)), 'equidistant', lonmin, lonmax, latmin, latmax, 7, 0, 16, 'biome',1)
map_plot_new(lon, lat, squeeze(field(6,:,:)), 'equidistant', lonmin, lonmax, latmin, latmax, 8, 300, 450, 'data_taka',1,2)
%map_plot(lon, lat, squeeze(field(10,:,:)), 'equidistant', lonmin, lonmax, latmin, latmax, 9, 0, 100, 'wind',1)
%map_plot(lon, lat, squeeze(field(11,:,:)), 'equidistant', lonmin, lonmax, latmin, latmax, 10, -2, 2, 'ssh',1)
%map_plot(lon, lat, squeeze(field(12,:,:)), 'equidistant', lonmin, lonmax, latmin, latmax, 11, 0, 1, 'curr speed',1)
map_plot_new(lon, lat, squeeze(field(12,:,:)), 'equidistant', lonmin, lonmax, latmin, latmax, 12, 280, 440, 'pCO2',1,2)
%map_plot(lon, lat, squeeze(field(14,:,:)), 'equidistant', lonmin, lonmax, latmin, latmax, 13, -20, 20, 'fCO2_anom',1)
map_plot_new(lon, lat, squeeze(field(7,:,:)), 'equidistant', lonmin, lonmax, latmin, latmax, 14, -2, 2, 'sst_anom',1,2)
map_plot_new(lon, lat, squeeze(field(8,:,:)), 'equidistant', lonmin, lonmax, latmin, latmax, 15, -2, 2, 'mld_anom',1,2)
map_plot_new(lon, lat, squeeze(field(9,:,:)), 'equidistant', lonmin, lonmax, latmin, latmax, 16, -2, 2, 'chl_anom',1,2)
map_plot_new(lon, lat, squeeze(field(10,:,:)), 'equidistant', lonmin, lonmax, latmin, latmax, 17, -2, 2, 'sss_anom',1,2)
map_plot_new(lon, lat, squeeze(field(11,:,:)), 'equidistant', lonmin, lonmax, latmin, latmax, 18, 0, 100, 'atm_anom',1,2)
map_plot_new(lon, lat, squeeze(field(13,:,:)), 'equidistant', lonmin, lonmax, latmin, latmax, 19, 0, 16, 'biomes',1,2)