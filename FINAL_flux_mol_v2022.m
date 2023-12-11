% This routine calculates the final flux estimates

%Peter Landschutzer 17.01.2012
%University of East Anglia, Norwich
%Edited by Annika Jersild July 2022

clear all
Input_Training_and_Labelling_new_GUI_v2022

%--------------------------------------------------------------------------
% 1) load all necessary data for calculation
%--------------------------------------------------------------------------

load('training_data/Ice/ice_v2022.mat');
load('training_data/wind/ERA5_component_wind_v2022.mat');

load('training_data/sst/sst_v2022.mat');
load('training_data/salinity/sss_v2022.mat');
load('training_data/pressure/pres_v2022.mat');
load('training_data/atm_co2/atm_pco2_grid_v2022.mat');
atm_co2=data_all;


% w=zeros(size(sst));
% w(:,:,:)=NaN;
% timevecXXXCC=1980:1/12:2014-1/12;
% ind=find(timevecXXXCC>=1990 & timevecXXXCC <2012);
% w(ind,:,:)=wind;
% clear wind
% wind=w;
% clear w

% [atm_co2]=calc_pco2_from_xco2(atm_co2,pressure,sst,sss);

% load the area

clear latsss lonsss timesss latwind lonwind sst_anom sss_anom

load('area.mat');

%--------------------------------------------------------------------------
% 2) calculate Schmidt number (Wanninkhof, 2014)
%--------------------------------------------------------------------------

A=2116.8;
B=136.25;
C=4.7353;
D=0.092307;
E=0.0007555;

Schmidt= A-B.*sst+C.*sst.^2-D.*sst.^3+E.*sst.^4;

clear A B C D E
%--------------------------------------------------------------------------
% 3) calculate solubility in mol/kg/atm (Weiss, 1974)
%--------------------------------------------------------------------------

A1=-58.0931;
A2=90.5069;
A3=22.2940;
B1=0.027766;
B2=-0.025888;
B3=0.0050578;

sst_abs=sst+273.15;

ln_sol=A1+A2*(100./sst_abs)+A3*log(sst_abs./100)+sss.*(B1+B2*(sst_abs./100)+B3*(sst_abs./100).^2);

sol=exp(ln_sol);

%mol/kg/atm is about mol/dm3/atm => calculate to mol/m3/mu atm
%try im mmol/m3/atm
%sol=sol.*10^6;
%in mmol/m3/mu atm
%sol=sol
% and in mol/m3/mu atm
sol=sol*10^-3;

clear sss ln_sol A1 A2 A3 B1 B2 B3 sst_abs latsst lonsst timesst


%--------------------------------------------------------------------------
% 5) load estimates and taka and calculate flux
%--------------------------------------------------------------------------

%pco2=nc_varget('/home/lapeter/Documents/post-doc/SOM-FFN/SOCOM-Analysis/SOCOM-fields/spco2_1986_2011_JENA_oc_v1.2.nc','spCO2');
% 
% data_all=zeros(408,180,360);
% 
% data_all(:,:,:)=NaN;
% data_all(61:396,:,:)=pco2;
% clear pco2;


if(strcmp(FFN_go,'yes')==1)
cd_dir_get=['output/NNoutput_SOCAT/' layer2take '/pCO2_' net2take '_' num2str(nnnumber) '.mat'];
elseif(strcmp(CL_go,'yes')==1)
cd_dir_get=['output/CLoutput_SOCAT/pCO2_' clusternr '.mat'];
elseif(strcmp(SOM_go,'yes')==1)
cd_dir_get=['output/SOMoutput_SOCAT/pCO2_' SOMnr '_' num2str(maplength) 'x' num2str(maphight) '.mat'];
elseif(strcmp(MLR_go,'yes')==1)
cd_dir_get=['output/MLRoutput_SOCAT/pCO2_' MLRnr '.mat'];
elseif(strcmp(BIOME_go,'yes')==1)
cd_dir_get=['output/BIOMEoutput_SOCAT/pCO2_' net2take '_' num2str(nnnumber) '.mat'];
elseif(strcmp(SSOM_go,'yes')==1)
cd_dir_get=['output/SUPERSOMoutput_SOCAT/' SSOMnr '_' num2str(nnnumber) '.mat'];
end
load(cd_dir_get);

% load('//ueahome1/env/hxh11bzu/.PC.USER.files/NTProfile/Documents/SOM/my_SOM/training_data/taka/Taka_pCO2_eth.mat')
% data_all=data_taka;
oc_co2_BIOME=data_all;
%oc_high=data_high;
%oc_low=data_low;
%oc_un=oc_high-oc_co2_BIOME;
[oc_co2_BIOME]=calc_pco2_from_fco2(oc_co2_BIOME,sst,pressure);
clear sst pressure
dfco2_BIOME=oc_co2_BIOME-atm_co2;
% end
%dfco2_high=oc_high-atm_co2;
%dfco2_low=oc_low-atm_co2;

% saveloc=['output/dpCO2_' net2take '_biome_' num2str(nnnumber) '_GV.mat'];
% save(saveloc, 'dfco2_BIOME');%, 'flux_est_low', 'flux_est_high');

%--------------------------------------------------------------------------
% 4) calculate transfer velocity in cm/hr
%--------------------------------------------------------------------------

transf_veloc_sw=0.251.*wind.*(Schmidt./660).^(-0.5);
%K = .251 based on Wannikhof 2014

load('area.mat')

transf_veloc_sw(isnan(oc_co2_BIOME))=NaN;
x_sw=16.5/(nanmean(area_average(transf_veloc_sw,area)));
transf_veloc_sw=transf_veloc_sw.*x_sw;
%and in m/yr:

transf_veloc_sw_new=transf_veloc_sw*24*365/100;

clear wind Schmidt


%--------------------------------------------------------------------------
% 6) calculate the flux in mol/m2/yr
%--------------------------------------------------------------------------

flux_est= transf_veloc_sw_new.*sol.*(1-seaice).*dfco2_BIOME;
%flux_est= transf_veloc_sw.*sol.*(1-seaice).*dfco2_BIOME;

save('output/parameter/parameterv2022.mat', 'sol', 'seaice', 'dfco2_BIOME', 'atm_co2', 'transf_veloc_sw');
%flux_est_low=transf_veloc_sw_new.*sol.*(1-seaice).*dfco2_low;
%flux_est_high=transf_veloc_sw_new.*sol.*(1-seaice).*dfco2_high;

% k_un=(((transf_veloc_sw_new).*0.3)./transf_veloc_sw_new).^2;
% oc_uncertainty=(oc_un./dfco2_BIOME).^2;
% 
% flux_un=sqrt(k_un+oc_uncertainty).*abs(flux_est);
% %flux_un=sqrt(oc_uncertainty).*abs(flux_est);

%[flux_estxx]=crop_field(flux_est,invareamin,invareamax);
%[flux_estxx]=crop_lat(flux_estxx,latcropmin,latcropmax);
%[flux_estxx]=crop_lon(flux_estxx,loncropmin,loncropmax);

%[flux_estxxx]=crop_field(oc_co2_BIOME,invareamin,invareamax);
%[flux_estxxx]=crop_lat(flux_estxxx,latcropmin,latcropmax);
%[flux_estxxx]=crop_lon(flux_estxxx,loncropmin,loncropmax);

% [flux_estxxxx]=crop_field(flux_un,invareamin,invareamax);
% [flux_estxxxx]=crop_lat(flux_estxxxx,latcropmin,latcropmax);
% [flux_estxxxx]=crop_lon(flux_estxxxx,loncropmin,loncropmax);

flux_estxx=flux_est;
flux_estxxx=oc_co2_BIOME;
%flux_estxxxx=flux_un;

test2=squeeze(nanmean(flux_estxx));
test3=squeeze(nanmean(flux_estxxx));
% test4=squeeze(nanmean(flux_estxxxx));

lat =-89.5:1:89.5;
lon =-179.5:1:179.5;
xx=size(lat);
yy=size(lon);
for ii=1:1:yy(2)
    
   latt(:,ii) = lat(:);   
end
for ii=1:1:xx(2)
   lonn(ii,:) = lon(:);  
end

%figure(1);
%contour_plot(lonn, latt, test2, 'equidistant', -180, 180, -90, 90, -5, 5, '')

map_plot_new(lon, lat, test2, 'equidistant', -180, 180, -90, 90, 1, -5, 5, 'flux',2,2)


% map_plot_new(lon, lat, squeeze(nanmean(dfco2_BIOME)), 'equidistant', -180, 180, -90, 90, 2, -50, 50, 'delta pCO2',2,2)
% printfilename=['dpCO2.png'];
% set(gcf,'PaperUnits','inches','PaperPosition',[0 0 16 7])
% print('-dpng', printfilename ,'-r200')
%--------------------------------------------------------------------------
% 7) save output
%--------------------------------------------------------------------------

saveloc=['output/flux/flux_' net2take '_biome_' num2str(nnnumber) '.mat'];
save(saveloc, 'flux_est');%, 'flux_est_low', 'flux_est_high');

%saveloc=['output/delta/dpCO2_' net2take '_biome_' num2str(nnnumber) '_GV.mat'];
%save(saveloc, 'dfco2_BIOME');%, 'flux_est_low', 'flux_est_high');

% saveloc=['output/flux/flux_taka.mat'];
% save(saveloc, 'flux_est');

