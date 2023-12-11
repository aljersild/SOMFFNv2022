% SMOOTH data

%Peter Landschutzer 17.01.2012
%University of East Anglia, Norwich
%Edited by Annika Jersild July 2022

clear all
Input_Training_and_Labelling_new_GUI_v2022

%--------------------------------------------------------------------------
% 1) load estimates
%--------------------------------------------------------------------------

if(strcmp(BIOME_go,'yes')==1)
cd_dir_get=['output/BIOMEoutput_SOCAT/pCO2_' net2take '_biome_' num2str(nnnumber) '.mat'];
elseif(strcmp(SSOM_go,'yes')==1)
cd_dir_get=['output/SUPERSOMoutput_SOCAT/' SSOMnr '_' num2str(nnnumber) '.mat'];
end

%load('output/BIOMEoutput_SOCAT/fco2_30yr.mat');

%temp1=fco2;
load(cd_dir_get);
temp1=data_all;


temp2=temp1;

ff=size(temp1);
cc=squeeze(nanmean(temp1));

new=sffilt('nanmean',temp2,[3,3,3]);


for rr=1:ff(1)
    newx=squeeze(new(rr,:,:));
    newx(isnan(cc))=NaN;
    newx(newx==0)=NaN;
    
    new2(rr,:,:)=newx;
    clear newx
    
end

clear new
new=new2;
clear new2

test1=squeeze(nanmean(new(25:end,:,:)));
test2=squeeze(nanmean(temp2(25:end,:,:)));
%diff=test1-test2;
lat1 =-89.5:1:89.5;
lon1 =-179.5:1:179.5;
[lon lat]=meshgrid(lon1,lat1);
time=min(year_output):1/12:max(year_output)+11/12;

map_plot_new(lon1, lat1, test1, 'equidistant', lonmin, lonmax, latmin, latmax,1, 280, 440, 'pCO2',1,2)
% printfilename=['pCO2_SOCCOM_2014_2015_only'];
% set(gcf,'PaperUnits','inches','PaperPosition',[0 0 16 7])
% print('-dpng', printfilename ,'-r200')


map_plot_new(lon1, lat1, test2, 'equidistant', lonmin, lonmax, latmin, latmax,2, 280, 440, 'pCO2',1,2)
%map_plot(lon, lat, diff, 'equidistant', lonmin, lonmax, latmin, latmax,3, -2, 2, 'pCO2 diff',1)



clear data_all
data_all=new;

%end
saveloc=['output/BIOMEoutput_SOCAT/pCO2_' net2take '_smoothed_biome_' num2str(nnnumber) '.mat'];
save(saveloc, 'data_all','data_biomes');

%save('output/BIOMEoutput_SOCAT/fco2_smoothed.mat', 'data_all');

