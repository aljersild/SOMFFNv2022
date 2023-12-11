%This routine can be used to combine 2 runs e.g. one with and one without
%chl to create a combine field

%Peter Landschutzer 17.01.2012
%University of East Anglia, Norwich
%edited by Annika Jersild July 2022

clear all
Input_Training_and_Labelling_new_GUI_v2022
%--------------------------------------------------------------------------
% 1) load estimates run1
%--------------------------------------------------------------------------

%run1
if(strcmp(FFN_go,'yes')==1)
cd_dir_get=['output/NNoutput_SOCAT/' layer2take '/pCO2_' net2take '_' num2str(nnnumber) '.mat'];
elseif(strcmp(CL_go,'yes')==1)
cd_dir_get=['output/CLoutput_SOCAT/pCO2_' clusternr '.mat'];
elseif(strcmp(SOM_go,'yes')==1)
cd_dir_get=['output/SOMoutput_SOCAT/pCO2_' SOMnr '_' num2str(maplength) 'x' num2str(maphight) '.mat'];
elseif(strcmp(MLR_go,'yes')==1)
cd_dir_get=['output/MLRoutput_SOCAT/pCO2_' MLRnr '.mat'];
elseif(strcmp(BIOME_go,'yes')==1)
cd_dir_get=['output/BIOMEoutput_SOCAT/pCO2_' net2take '_biome_' num2str(nnnumber) '.mat'];
elseif(strcmp(SSOM_go,'yes')==1)
cd_dir_get=['output/SUPERSOMoutput_SOCAT/' SSOMnr '_' num2str(nnnumber) '.mat'];
end
%run2

%--------------------------------------------------------------------------
% 2) load estimates run2
%--------------------------------------------------------------------------

prompt={'net2compare','nnnumber1'};
if(strcmp(FFN_go,'yes')==1)
defans1={net2take,'nnnumber'};
elseif(strcmp(CL_go,'yes')==1)
defans1={clusternr,'nnnumber'};
elseif(strcmp(SOM_go,'yes')==1)
defans1={SOMnr,'nnnumber'};
elseif(strcmp(MLR_go,'yes')==1)
defans1={MLRnr,'nnnumber'};
elseif(strcmp(BIOME_go,'yes')==1)
defans1={net2take,'nnnumber'};
elseif(strcmp(SSOM_go,'yes')==1)
defans1={SSOMnr,'nnnumber'};
end
fields = {'net2compare','nnnumber1'}; 
options.Resize='on';
options.WindowStyle='modal';
options.Interpreter='none';
info = inputdlg(prompt, 'compare with', 1, defans1,options);
if ~isempty(info)              %see if user hit cancel
   info = cell2struct(info,fields);
   net2compare    = info.net2compare;
   nnnumber1    = info.nnnumber1;
end      
clear prompt defans fields info

if(strcmp(FFN_go,'yes')==1)
cd_dir_get1=['output/NNoutput_SOCAT/' layer2take '/pCO2_' net2compare '_' num2str(nnnumber1) '.mat'];
elseif(strcmp(CL_go,'yes')==1)
cd_dir_get1=['output/CLoutput_SOCAT/pCO2_' net2compare '.mat'];
elseif(strcmp(SOM_go,'yes')==1)
cd_dir_get1=['output/SOMoutput_SOCAT/pCO2_' net2compare '_' num2str(maplength) 'x' num2str(maphight) '.mat'];
elseif(strcmp(MLR_go,'yes')==1)
cd_dir_get1=['output/MLRoutput_SOCAT/pCO2_' net2compare '.mat'];
elseif(strcmp(BIOME_go,'yes')==1)
cd_dir_get1=['output/BIOMEoutput_SOCAT/pCO2_' net2compare '_biome_' num2str(nnnumber1) '.mat'];
elseif(strcmp(SSOM_go,'yes')==1)
cd_dir_get1=['output/SUPERSOMoutput_SOCAT/' net2compare '_' num2str(nnnumber1) '.mat'];
end

load(cd_dir_get);
% creating a temp file just in case sth goes wrong.
sL = ['data_all'];
temp=eval(sL);
sT = ['data_biomes'];
temp20=eval(sT);

load(cd_dir_get1);
% creating a temp file just in case sth goes wrong.
sL1 = ['data_all'];
temp1=eval(sL1);
sT1 = ['data_biomes'];
temp21=eval(sT1);

clear data_all data_biomes
%--------------------------------------------------------------------------
% 3) combine and save
%--------------------------------------------------------------------------

fco2=temp;
ind=find(isnan(fco2));
fco2(ind)=temp1(ind);
clear ind
biomes=temp20;
ind=find(isnan(biomes));
biomes(ind)=temp21(ind);
%biomes(ind)=temp21(ind)+1000;
%biomes(biomes>100)=17;

data_all=fco2;
data_biomes=biomes;

data_all(data_all<0)=NaN;

%--------------------------------------------------------------------------
% exclude Arctic
%--------------------------------------------------------------------------

lat=-89.5:1:89.5;
lon=-179.5:1:179.5;

[lonn latt]=meshgrid(lon,lat);

%crop everything >80N

% ind=find(latt>80);
% data_all(:,ind)=NaN;
% data_biomes(:,ind)=NaN;
% 
% 
% %crop Arctic everywhere else
% 
% ind2=find(lonn>30 & latt>65);
% ind3=find(lonn<-90 & latt>65);
% data_all(:,ind2)=NaN;
% data_biomes(:,ind2)=NaN;
% data_all(:,ind3)=NaN;
% data_biomes(:,ind3)=NaN;

fco2=squeeze(nanmean(data_all(end-11:end,:,:)));
fco21=squeeze(nanmean(data_all(end-23:end-12,:,:)));

% fco2=squeeze(nanmean(data_all(120:131,:,:)));
% fco21=squeeze(nanmean(data_all(132:143,:,:)));

map_plot_new(lon, lat, fco2, 'equidistant', lonmin, lonmax, latmin, latmax,2, 280, 600, '',1,1)
map_plot_new(lon, lat, fco21, 'equidistant', lonmin, lonmax, latmin, latmax,3, 280, 450, '',1,1)


%--------------------------------------------------------------------------
% save
%--------------------------------------------------------------------------

saveloc=['output/BIOMEoutput_SOCAT/pCO2_' net2take '_' net2compare '-combined_biome_' num2str(nnnumber) '.mat'];
save(saveloc, 'data_all','data_biomes');

%clear fco2 temp temp1 temp2 temp3 temp4 temp5 ind ind2 ind3 ind1 ind4 sL5 sL sL1 sL2 sL3 sL4
%clear all
