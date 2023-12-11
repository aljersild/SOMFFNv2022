% Welcome to the simple things: this code uses the NN matlab package and
% creates monthly generated basin wide fCO2 maps derrived from the SOM. The
% input file specifies the year and month and the lat-lon area. 

%Peter Landschutzer 17.01.2012
%University of East Anglia, Norwich
%Edited by Annika Jersild July 2022
clear all
Input_Training_and_Labelling_new_GUI_v2022

%--------------------------------------------------------------------------
% 1) load cluster data 
%--------------------------------------------------------------------------

load('training_data/mld/mld_clim_v2022.mat');
clear lonmld latmld
load('training_data/salinity/sss_v2022.mat');
clear lonsss latsss
load('training_data/taka/Taka_pCO2_eth_v2022.mat');
load('training_data/sst/sst_v2022.mat');

%--------------------------------------------------------------------------
% 2) take 20-year average
%--------------------------------------------------------------------------

timevec=1980:1/12:2022-1/12;
ind=find(timevec>=1982 &timevec<2022);

sst=sst(ind,:,:);
mld=mld(ind,:,:);
sss=sss(ind,:,:);
data_taka=data_taka(ind,:,:);


for uu=1:12
   sst_an(uu,:,:)=nanmean(sst(uu:12:end,:,:));
   mld_an(uu,:,:)=nanmean(mld(uu:12:end,:,:));
   sss_an(uu,:,:)=nanmean(sss(uu:12:end,:,:));
   taka_an(uu,:,:)=nanmean(data_taka(uu:12:end,:,:));
end

lonsst_an=lonsst(1:12,:,:);
latsst_an=latsst(1:12,:,:);
monthx=zeros(size(latsst_an));

for ww=1:12
    monthx(ww,:,:)=ww;
end

%--------------------------------------------------------------------------
% 3) reshape and rearrange for SOM
%--------------------------------------------------------------------------


sst_new = reshape(sst_an, prod(size(sst_an)), 1);
sss_new = reshape(sss_an, prod(size(sss_an)), 1);
mld_new = reshape(mld_an, prod(size(mld_an)), 1);
taka_new = reshape(taka_an, prod(size(taka_an)), 1);
lat_new = reshape(latsst_an, prod(size(latsst_an)), 1);
lon_new = reshape(lonsst_an, prod(size(lonsst_an)), 1);
month_new = reshape(monthx, prod(size(monthx)), 1);

clear lon lat month

temp2=[sst_new mld_new sss_new taka_new lat_new lon_new month_new];


[index2,Train_SOM]=nanremoveV3_SOM(temp2);

LAT2=lat_new(index2);
LON2=lon_new(index2);
month=month_new(index2);
%Label_SOM=fCO2;

clear index1 index2 temp1 temp2 sum_index sum_index1


%--------------------------------------------------------------------------
% 4) normalize data 
%--------------------------------------------------------------------------

%net = train(net,pn,tn);
%log normalization
%[Label_SOM, Train_SOM]=log_norm(Label_SOM,Train_SOM,idx1);

% Normalization using som toolbox not idx adjusted!!!
%[Label_SOM, Train_SOM, Label_FFN, Train_FFN]=norm_Maciej_SF(Label_SOM,Train_SOM,Label_FFN,Train_FFN);

%normalization using mapminmax
%[Label_SOM, Train_SOM]=norm_minmax(Label_SOM,Train_SOM);

% if(~isempty(find(idx1==14, 1)));
% inddxx=find(idx1==14);
% Label_SOM(:,inddxx)=Label_SOM(:,inddxx)*3;
% Train_SOM(:,inddxx)=Train_SOM(:,inddxx)*3;    
% clear inddxx
% end

clear idx1 idx2
%--------------------------------------------------------------------------
% 5) SOM part to identify biomes
%--------------------------------------------------------------------------

    
net = selforgmap([maplength maphight],100,3,'hextop','dist');
net.trainParam.epochs = epochnr;
[net,tr] = train(net,Train_SOM');
    
y = net(Train_SOM');
classes = vec2ind(y);    

%--------------------------------------------------------------------------
% 6) Smoothing of biomes
%--------------------------------------------------------------------------


for month2plot=1:12
 
            index4=find(month==month2plot);

            classesYYY=classes(index4);

                lat =-89.5:1:89.5;
                lon =-179.5:1:179.5;
                LATY = LAT2(index4);
                LONY = LON2(index4);
                [classesY]=vec_to_array2(classesYYY,lon,lat,LONY,LATY);

                biomeY=classesY;
                biomeY2(month2plot,:,:)=biomeY(:,:);
end

lon=-179.5:1:179.5;
lat=-89.5:1:89.5;


%--------------------------------------------------------------------------
% 5.5) arctic 1 biome
%--------------------------------------------------------------------------

[lonn latt]=meshgrid(lon,lat);
%crop Arctic everywhere else

ind2=find(lonn>30 & latt>65);
ind3=find(lonn<-90 & latt>65);
ind4=find(latt>80);

indtot2=[ind2; ind3; ind4];
indtot=unique(indtot2);

biomeset=biomeY2(:,indtot);
biomeset(isnan(biomeset))=[];
mode_b=mode(biomeset);


biomeY2(:,indtot)=mode_b;

%--------------------------------------------------------------------------
% no smoothing
%--------------------------------------------------------------------------



%map_plot(lon, lat, squeeze(biomeY2(10,:,:)), 'equidistant', -180, 180, -90, 90,2, 0, 16, 'pCO2',1,1)

gg2=size(biomeY2)
         
    for yy=1:gg2(1)
        yy
        new=sffilt('mode',squeeze(biomeY2(yy,:,:)),[10 10]);
        for bio=1:16
            ind=find(squeeze(biomeY2(yy,:,:))~=bio);
            bio_new=squeeze(biomeY2(yy,:,:));
            bio_new(ind)=0;
            bio_new(bio_new==bio)=1;
            BW2 = bwareaopen(bio_new, 15);
            BW2=double(BW2);
            BW2(BW2==0)=NaN;
            BW2(~isnan(BW2))=bio;
            if bio==1
                BW=BW2;
            else
                ind2= isnan(BW);
                BW(ind2)=BW2(ind2);
            end
            
        clear BW2 bio_new ind ind2
        end
     new2(yy,:,:)=reshape(new,180,360);
     BX(yy,:,:)=BW;
     clear BW
    end  
    
    clear new
    new=new2;
    clear new2
    new(isnan(biomeY2))=NaN;    
    
    % map_plot(lon, lat, squeeze(new(10,:,:)), 'equidistant', -180, 180, -90, 90,4, 0, 16, 'pCO2',1,1)
    
    
    indg=find(isnan(BX));
    BX(indg)=new(indg);
    
    clear new
    new=BX;
    newL=new;


for ii=1:12:504
   biomes(ii:ii+11,:,:)=new(:,:,:); 
end

%--------------------------------------------------------------------------
% 7) save and plot 3-D biomes
%--------------------------------------------------------------------------


%array_test=squeeze(biomes(2,:,:));
array_test=squeeze(mode(biomes));

saveloc=['output/BIOMEoutput_SOCAT/networks/' SOMnr '_biome_' num2str(maplength) 'x' num2str(maphight) '.mat'];
save(saveloc, 'net', 'tr','biomes');


figure(6);clf; 
m_proj('equidistant','lon',[-180 180],'lat',[-90 90]);
h=m_pcolor(lon,lat,array_test);
set(h,'edgecolor','none');
colorbar('hori');caxis([0 16]);
%colormap(redblue);
m_coast('patch','w');
m_grid('xtick',[ -135 -90 -45 0 45 90 135], ...
    'ytick',[-90 -60 -30 0 30 60 90],'fontsize',10);



