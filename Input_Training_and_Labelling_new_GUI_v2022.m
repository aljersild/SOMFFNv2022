% This is new: a flexible input file where LON LAT YR and MONTH for plots
% are fixed. All related codes of the NN use this input file, so the
% results are consistent. some other parameters for the NN training are
% adjusted in this file as well. And the best: its all GUI!!!

%Peter Landschutzer 17.01.2012
%University of East Anglia, Norwich

%====================== add path to some tools ============================
addpath /Users/annikajersild/Documents/MATLAB/m_map   % mapping tool for plots
addpath /Users/annikajersild/Documents/MATLAB/mexcdf/mexnc % to read file info and variables in matlab
addpath /Users/annikajersild/Documents/MATLAB/mexcdf/snctools % idem dito
addpath /Users/annikajersild/Documents/MATLAB/nansuite % mean ignoring nan
%addpath /Users/annikajersild/Documents/MATLAB/functions
addpath /Users/annikajersild/Documents/Research_2022/UncertaintyAnalysis/functions
addpath /Users/annikajersild/Documents/MATLAB
%============== control parameter =========================================

prompt={'year_min','year_max', 'year2find', 'month2plot'};
defans={'1982','2020', '2005', '1'};
fields = {'year_min','year_max', 'year2find', 'month2plot'};
%options.Resize='on';
%options.WindowStyle='modal';
%options.Interpreter='none';
%info='a';
%info = inputdlgcol(prompt, 'please enter time and geographic boarders', 1, defans,options,2);
info = inputdlg(prompt, 'please enter time', 1, defans);
if ~isempty(info)              %see if user hit cancel
   info = cell2struct(info,fields);
   year_min   = str2num(info.year_min);
   year_max   = str2num(info.year_max);
   year2find  = str2num(info.year2find);
   month2plot = str2num(info.month2plot);
   %xx=str2num(info.xx);
end
clear prompt defans fields info
year_output=(year_min:1:year_max);
timevec=1982:1/12:year_max+1;


prompt={'latmin', 'latmax', 'lonmin', 'lonmax','invareamin','invareamax', 'latplotmin', 'latplotmax', 'lonplotmin', 'lonplotmax', 'version'};
defans={'-90', '90', '-180', '180','1','11','-90', '90', '-180', '180', '2'};
fields = {'latmin', 'latmax', 'lonmin', 'lonmax','invareamin','invareamax', 'latplotmin', 'latplotmax', 'lonplotmin', 'lonplotmax', 'version'};
%info='a';
%info = inputdlgcol(prompt, 'please enter time and geographic boarders', 1, defans,options,2);
info = inputdlg(prompt, 'please enter geographic boarders', 1, defans);
if ~isempty(info)              %see if user hit cancel
   info = cell2struct(info,fields);
   loncropmin = str2num(info.lonmin);%convert string to number
   loncropmax = str2num(info.lonmax);
   latcropmin = str2num(info.latmin);   
   latcropmax = str2num(info.latmax);
   invareamin = str2num(info.invareamin);
   invareamax = str2num(info.invareamax);
   lonmin     = str2num(info.lonplotmin);%convert string to number
   lonmax     = str2num(info.lonplotmax);
   latmin     = str2num(info.latplotmin);   
   latmax     = str2num(info.latplotmax);
   version    = str2num(info.version);
   %xx=str2num(info.xx);
end
clear prompt defans fields info

%lat
%lon
%time

prompt={'FFN','SOM', 'Clustering', 'MLR', 'BIOME', 'SSOM'};
defans={'no','no', 'no', 'no', 'yes', 'no'};
fields = {'FFN','SOM', 'Clustering', 'MLR', 'BIOME', 'SSOM'};
info = inputdlg(prompt, 'please enter method', 1, defans);
if ~isempty(info)              %see if user hit cancel
   info = cell2struct(info,fields);
   FFN_go   = info.FFN;
   SOM_go   = info.SOM;
   CL_go    = info.Clustering;%convert string to number
   MLR_go   = info.MLR;%convert string to number
   BIOME_go = info.BIOME;
   SSOM_go  = info.SSOM;
end
clear prompt defans fields info

if(strcmp(FFN_go,'yes')==1)||(strcmp(BIOME_go,'yes')==1)
prompt={'net2take','netlayer', 'layer2take', 'nnnumber','load_trained_net'};
defans={'GCBv2020_GCBv2020a-combined_smoothed','2layer', '2layer_plus_coast_v1p5', '50','"yes or no"'};
fields = {'net2take','netlayer', 'layer2take', 'nnnumber','load_trained_net'};
info = inputdlg(prompt, 'please enter method', 1, defans);
if ~isempty(info)              %see if user hit cancel
   info = cell2struct(info,fields);
   net2take         = info.net2take;
   netlayer         = info.netlayer;
   layer2take       = info.layer2take;%convert string to number
   nnnumber         = str2num(info.nnnumber);%convert string to number
   load_trained_net = info.load_trained_net;
end
clear prompt defans fields info
end

if(strcmp(SOM_go,'yes')==1)||(strcmp(BIOME_go,'yes')==1)
prompt={'SOMnr','maplength', 'maphight', 'epochnr','load_trained_SOM','hc clusters'};
defans={'SOM','7', '7', '5','"yes or no"', '10'};
fields = {'SOMnr','maplength', 'maphight', 'epochnr','load_trained_SOM','hc_clusters'};
info = inputdlg(prompt, 'please enter method', 1, defans);
if ~isempty(info)              %see if user hit cancel
   info = cell2struct(info,fields);
   SOMnr            = info.SOMnr;
   maplength        = str2num(info.maplength);
   maphight         = str2num(info.maphight);%convert string to number
   epochnr          = str2num(info.epochnr);%convert string to number
   load_trained_SOM = info.load_trained_SOM;
   hc_clusters      = str2num(info.hc_clusters);
end      
clear prompt defans fields info
end

if(strcmp(CL_go,'yes')==1)
prompt={'clusternr','nb_cluster'};    
defans={'cl100','100'};
fields = {'clusternr','nb_cluster'};  
info = inputdlg(prompt, 'please enter method', 1, defans);
if ~isempty(info)              %see if user hit cancel
   info = cell2struct(info,fields);
   clusternr    = info.clusternr;
   nb_cluster   = str2num(info.nb_cluster);
end      
clear prompt defans fields info
end

if(strcmp(MLR_go,'yes')==1)
prompt={'MLRnr'};    
defans={'MLR100'};
fields = {'MLRnr'};  
info = inputdlg(prompt, 'please enter method', 1, defans);
if ~isempty(info)              %see if user hit cancel
   info = cell2struct(info,fields);
   MLRnr    = info.MLRnr;
end      
clear prompt defans fields info
end

if(strcmp(SSOM_go,'yes')==1)
prompt={'SSOMnr', 'nnnumber'};    
defans={'fCO2_Super_SOM_biome_1', '2'};
fields = {'SSOMnr', 'nnnumber'};  
info = inputdlg(prompt, 'please enter method', 1, defans);
if ~isempty(info)              %see if user hit cancel
   info = cell2struct(info,fields);
   SSOMnr    = info.SSOMnr;
   nnnumber  = info.nnnumber;
end      
clear prompt defans fields info
end




