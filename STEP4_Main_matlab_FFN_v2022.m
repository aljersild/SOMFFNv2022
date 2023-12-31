% Welcome to the simple things: this code uses the NN matlab package and
% creates monthly generated basin wide fCO2 maps derrived from the SOM. The
% input file specifies the year and month and the lat-lon area. 

%Peter Landschutzer 17.01.2012
%University of East Anglia, Norwich
%Edited by Annika Jersild July 2022


clear all
Input_Training_and_Labelling_new_GUI_v2022
%fileID = fopen(['protocol_' net2take '.txt'],'w');
%--------------------------------------------------------------------------
% 1) open GUI to define basin
%--------------------------------------------------------------------------

prompt={'ocean','deactivate_testnum','val_check_only','performance_function'};    
defans={'please enter "Atlantic" or "Global"','"yes" or "no"','yes or no','mse or sse'};
fields = {'ocean','deactivate_testnum','val_check_only','performance_function'};  
info = inputdlg(prompt, 'please enter if global or Atlantic', 1, defans);
if ~isempty(info)              
   info = cell2struct(info,fields);
   ocean  = info.ocean;
   deactivate_testnum  = info.deactivate_testnum;
   val_check_only  = info.val_check_only;
   performance_function = info.performance_function;
end      
clear prompt defans fields info

%--------------------------------------------------------------------------
% 2) load the labelling and training data
%--------------------------------------------------------------------------

temp3=0;

for year=min(year_output):max(year_output)
%for SOCAT
if(strcmp(ocean,'Global')==1)
path1 =(['output/Ldata/Ldata_' num2str(year) '_v2022.mat']);
elseif(strcmp(ocean,'Atlantic')==1)
path1 =(['output/Ldata/Ldata_' num2str(year) '_coincided_v2.mat']);
end
load(path1);
% creating a temp file just in case sth goes wrong.
sL = ['Ldata_' num2str(year) '_coincided_v2'];
temp1=eval(sL);
[temp1]=clear_datasetV2(temp1, 'yes', latcropmin, latcropmax,'yes', loncropmin, loncropmax, 'no', 'no', 'no');

if(year==min(year_output)) %2002 for UEA and BATS
%if(year==2014)
    temp3=temp1;
else
    temp3=[temp3; temp1];
end

clear Ldata_1982_coincided_v2 Ldata_1983_coincided_v2 Ldata_2012_coincided_v2
clear Ldata_1984_coincided_v2 Ldata_1985_coincided_v2 Ldata_1986_coincided_v2
clear Ldata_1987_coincided_v2 Ldata_1988_coincided_v2  Ldata_1989_coincided_v2
clear path1 temp1 sL Ldata_1998_coincided_v2 Ldata_1999_coincided_v2
clear Ldata_1992_coincided_v2 Ldata_1993_coincided_v2 Ldata_1994_coincided_v2
clear Ldata_1995_coincided_v2 Ldata_1996_coincided_v2  Ldata_1997_coincided_v2
clear Ldata_2000_coincided_v2 Ldata_2001_coincided_v2 Ldata_2002_coincided_v2
clear Ldata_2003_coincided_v2 Ldata_2004_coincided_v2  Ldata_2005_coincided_v2
clear Ldata_2006_coincided_v2 Ldata_2007_coincided_v2 Ldata_2008_coincided_v2
clear Ldata_2009_coincided_v2 Ldata_2010_coincided_v2 Ldata_2011_coincided_v2
clear Ldata_2012_coincided_v2 Ldata_2013_coincided_v2 Ldata_2014_coincided_v2
clear Ldata_1990_coincided_v2 Ldata_1991_coincided_v2 Ldata_2015_coincided_v2
clear Ldata_2016_coincided_v2 Ldata_2017_coincided_v2 Ldata_2018_coincided_v2
clear Ldata_2019_coincided_v2
end
temp1=temp3;
clear temp3;
for year=min(year_output):max(year_output)
if(strcmp(ocean,'Global')==1)
path2 =(['output/Tdata/Tdata_' num2str(year) '_v2022.mat']);
elseif(strcmp(ocean,'Atlantic')==1)
path2 =(['output/Tdata/Tdata_' num2str(year) 'V2.mat']);
end
load(path2);
% creating a temp file just in case sth goes wrong.
sT = ['Tdata_' num2str(year)];
temp2=eval(sT);
[temp2]=clear_datasetV2(temp2, 'yes', latcropmin, latcropmax,'yes', loncropmin, loncropmax, 'no', 'no', 'no');

if(year==min(year_output)) %2002 for UEA and BATS
    temp4=temp2;
else
    temp4=[temp4; temp2];
end
clear path2 temp2 sT Tdata_1998 Tdata_1999 Tdata_2000 Tdata_2001 Tdata_2002 Tdata_2010 
clear Tdata_2003  Tdata_2004 Tdata_2005 Tdata_2006 Tdata_2007 Tdata_2008 Tdata_2009 Tdata_2011
clear Tdata_1992  Tdata_1993 Tdata_1994 Tdata_1995 Tdata_1996 Tdata_1997 Tdata_2012 Tdata_2013
clear Tdata_1982  Tdata_1983 Tdata_1984 Tdata_1985 Tdata_1986 Tdata_1987 Tdata_1988 Tdata_1989
clear Tdata_1990  Tdata_1991 Tdata_2015 Tdata_2016 Tdata_2017 Tdata_2018 Tdata_2019 Tdata_2014
end
temp2=temp4;
clear temp4;

%--------------------------------------------------------------------------
% 3) NAN's of both labelling and training data have to be removed
%--------------------------------------------------------------------------

[index1,index2,Label_FFN,Train_FFN,idx2]=nanremoveV3(temp1,temp2);

fCO2 = temp1(index1,16);       
Month1=temp1(index1,2);
Month2=temp2(index2,2);
Year2=temp2(index2,1);
Year1=temp1(index1,1);
LAT1 = temp2(index2,3);
LON1 = temp2(index2,4);
LAT2 = temp1(index1,3);
LON2 = temp1(index1,4);
classes1 = temp2(index2,17);
classes = temp1(index1,17);

% ind=find(fCO2>800);
% 
% fCO2(ind)=[];
% Month1(ind)=[];
% Month2(ind)=[];
% Year2(ind)=[];
% Year1(ind)=[];
% LAT1(ind)=[];
% LON1(ind)=[];
% LAT2(ind)=[];
% LON2(ind)=[];
% classes1(ind)=[];
% classes(ind)=[];
% Label_FFN(ind,:)=[];
% Train_FFN(ind,:)=[];
% 
% ind2=find(fCO2<50);
% 
% fCO2(ind2)=[];
% Month1(ind2)=[];
% Month2(ind2)=[];
% Year2(ind2)=[];
% Year1(ind2)=[];
% LAT1(ind2)=[];
% LON1(ind2)=[];
% LAT2(ind2)=[];
% LON2(ind2)=[];
% classes1(ind2)=[];
% classes(ind2)=[];
% Label_FFN(ind2,:)=[];
% Train_FFN(ind2,:)=[];


content = unique(classes);

clear temp1 temp2 sum_index sum_index1 data_biomes lonbiomes latbiomes index1 index2

%---------------------------------------------------------------------
%     classes1(classes1==16)=15;
%     classes(classes==16)=15;
%     content(content==16)=[];
% fprintf(fileID,'%6s/n','**********************************');
% fprintf(fileID,'%6s/n',[' normalized SOM weigths ' num2str(length(idx1)) 'parameters x ' num2str(maplength*maphight) 'neurons']);
% fprintf(fileID,'%6.2f %6.2f %6.2f %6.2f %6.2f/n',b);
% fprintf(fileID,'%6s/n','**********************************');

%--------------------------------------------------------------------------
% end clustering
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% 4) Backprop part for every Neuron
%--------------------------------------------------------------------------

content4SOM=content;
content(isnan(content))=[];

hits_FFN=0;
hits_mean=0;
hits_empty=0;
    
%--------------------------------------------------------------------------
% start FFN part for every biome
%--------------------------------------------------------------------------
hitlist=zeros(length(content),1);
hitdata=zeros(length(content),1);
for biome1=1:length(content)
    
    biome=content(biome1)

    ind_FFN=find(classes==biome);
    ind_FFN1=find(classes1==biome);
    hitdata(biome1)=length(ind_FFN);
    
    classes_Train=classes1(ind_FFN1);
    Train1=Train_FFN(ind_FFN1,:);
    Label1=Label_FFN(ind_FFN,:);
    fCO2_FFN=fCO2(ind_FFN);
    Month_Train_FFN=Month2(ind_FFN1);
    Year_Train_FFN=Year2(ind_FFN1);
    LAT_Train = LAT1(ind_FFN1);
    LON_Train = LON1(ind_FFN1);


    if(strcmp(deactivate_testnum,'yes')==1)
        testnum=100000000;
        nnnumber=200;
    else
        testnum=100;
    end

    %----------------------------------------------------------------------
    % if less than 100 observations per biome perform a MLR instead of FFN
    %----------------------------------------------------------------------
    
    if(length(ind_FFN)<testnum)&&(length(ind_FFN)>1)
        hits_mean=hits_mean+1;
        hitlist(biome1)=1;
        predictor=[ones(size(Label1),1) Label1];
        [b1,bint] =regress(fCO2_FFN, predictor);
        trainer=[ones(length(Train1),1) Train1];
        SOMCO2output=[trainer]*b1;
        
        for year2go=min(year_output):max(year_output)

        yearyeah(year2go-(min(year_output)-1))=year2go;
            for month2plot=1:12
        
            index3=find(Year_Train_FFN==year2go & Month_Train_FFN==month2plot);
            new_fCO2=SOMCO2output(index3);
            classes_Train1=classes_Train(index3);
            Train2 = Train1(index3,:);
            TF=isempty(Train2);
                if(TF==0)
                lat =-89.5:1:89.5;
                lon =-179.5:1:179.5;
                LAT = LAT_Train(index3);
                LON = LON_Train(index3);
                [field]=vec_to_array2(new_fCO2,lon,lat,LON,LAT);
                [classes_field]=vec_to_array2(classes_Train1,lon,lat,LON,LAT);
                fco2=field;
                biome_regimes1=classes_field;
                else
                fco2(1:180,1:360)=NaN;
                biome_regimes1(1:180,1:360)=NaN;
                end
            temp_array(month2plot,:,:)=fco2(:,:);
            biome_regimes2(month2plot,:,:)=biome_regimes1(:,:);
            end
                if(year2go==min(year_output))
                temp_array1=temp_array;
                biome_regimes3=biome_regimes2;
                else
                temp_array1=[temp_array1; temp_array];
                biome_regimes3=[biome_regimes3;biome_regimes2];
                end
        end
        
    %----------------------------------------------------------------------
    % if enough observations per biome perform FFN
    %----------------------------------------------------------------------
            
    elseif(length(ind_FFN)>testnum)
        %------------------------------------------------------------------
        %first decide on the hidden layer
        %------------------------------------------------------------------
        a=size(Label1);
        b=a(1);
        c=a(2);
        %rule of thumb: size is bigger than output layer (=1*a(2)) and
        %smaller than input layer = a(1)
        %therefor:
        nnminnum=2;
        % weigths to training data shoult be more than 30: maxnum*a(2)+maxnum <=a(1)
        % + again ~ 10 % are labelling data so 0.9;
        nnmaxnum1=(0.9*a(1))/30;
        nnmaxnum=floor(sqrt(nnmaxnum1));
        if nnmaxnum<2
            nnmaxnum=2;
        elseif nnmaxnum>14
            nnmaxnum=14;
        end
        %------------------------------------------------------------------
        
        hits_FFN=hits_FFN+1;
        hitlist(biome1)=2;
        
        % start pre-training with variable neuron number
        for nnnumber11=nnminnum:nnmaxnum
        %for nnnumber11=2:14
        %nnnumber=nnnumber11*10;
        nnnumber=ceil((nnnumber11*nnnumber11)/2);
        %nnnumber=50;
        net = feedforwardnet(nnnumber);
        
        
        Label1(:,end+1)=fCO2_FFN;         
        ix=size(Label1,2);
        sorted=sortrows(Label1,ix);
        clear Label1 fCO2_FFN ix
        Label1=sorted(:,1:end-1);
        fCO2_FFN=sorted(:,end);
        clear sorted
        
%       Amary criterion
        ixx=1/sqrt(2*nnnumber);
        net = configure(net,Label1',fCO2_FFN');

        net.divideParam.trainRatio = 1-ixx;
        net.divideParam.valRatio = ixx;
        net.divideParam.testRatio = 0;
 
        clear ix ix1 ixx

        %train the network
        [net,tr] = train(net, Label1', fCO2_FFN');
        %Add in code to not pop up all the windows
        net.trainParam.showWindow = false;
        
        checkout(nnnumber11-1)=abs(tr.vperf(tr.best_epoch+1)-tr.perf(tr.best_epoch+1));
        checkout2(nnnumber11-1)=nnnumber;
        val_check(nnnumber11-1)=tr.vperf(tr.best_epoch+1);
        perf_check(nnnumber11-1)=tr.perf(tr.best_epoch+1);
        end
        
        checkout_test=checkout+val_check;
        if(strcmp(val_check_only,'yes')==1)
        ind_checkout=find(val_check==min(val_check));
        else
        ind_checkout=find(checkout_test==min(checkout_test));
        end
        nnnumberX=checkout2(ind_checkout);
        val_checkX=val_check(ind_checkout);
        perf_checkX=perf_check(ind_checkout);
        
        
        for mruns=1:10
        net = feedforwardnet(nnnumberX);
        net = configure(net,Label1',fCO2_FFN');
        %net.trainParam.epochs = 50;
        ixx=1/sqrt(2*nnnumberX);
        net.performFcn=performance_function;
        net.divideParam.trainRatio = 1-ixx;
        net.divideParam.valRatio = ixx;
        net.divideParam.testRatio = 0;

        %train the network
        [net,tr] = train(net, Label1', fCO2_FFN');
        %saveloc=['output/BIOMEoutput_SOCAT/networks/' SOMnr '_biome_' num2str(maplength) 'x' num2str(maphight) '.mat'];
        %save(saveloc, 'net', 'tr');
        %make several simulations in mruns loop and take the average
%         a=net.IW;
%         b=cell2mat(a);
%         fprintf(fileID,'%6s/n','**********************************');
%         fprintf(fileID,'%6s/n',[' FFN weigths biome' num2str(biome1) ' ' num2str(length(idx2)) 'parameters x ' num2str(nnnumberX) 'neurons']);
%         fprintf(fileID,'%8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f/n',b);
%         fprintf(fileID,'%6s/n','**********************************');
        new_fCO2_sim = sim(net, Train1');
%      new_fCO2X =sim(net,Label1');

            %if(mruns==1)
            %fCO2_sim=new_fCO2_sim;
            %else
            %fCO2_sim=fCO2_sim+new_fCO2_sim;
            %end
            new_fCO2_sim(new_fCO2_sim<0)=NaN;
            fCO2_sim2(mruns,:)=new_fCO2_sim;           
            
        end
       %new_fCO2_sim=fCO2_sim/mruns;
       new_fCO2_sim=squeeze((nanmean(fCO2_sim2,1)));
       clear fCO2_sim2
       fCO2_sim=new_fCO2_sim;
       clear new_fCO2_sim
       
        clear checkout checkout2 val_check perf_check val_ind perf_ind

       %saveloc=['output/BIOMEoutput_SOCAT/networks/' net2take '_biome_' num2str(biome) '.mat'];
       %save(saveloc, 'net', 'tr');
        
    %----------------------------------------------------------------------
    % array loop
    %----------------------------------------------------------------------       
       
        for year2go=min(year_output):max(year_output)     

        yearyeah(year2go-(min(year_output)-1))=year2go;    
            for month2plot=1:12

            index3=find(Year_Train_FFN==year2go & Month_Train_FFN==month2plot);
            classes_Train1=classes_Train(index3);
            Train2 = Train1(index3,:);

            TF=isempty(Train2);
                if(TF==0)
                LAT = LAT_Train(index3);
                LON = LON_Train(index3);
                lat =-89.5:1:89.5;
                lon =-179.5:1:179.5;
                new_fCO2=fCO2_sim(index3);
                [field]=vec_to_array2(new_fCO2,lon,lat,LON,LAT);
                [classes_field]=vec_to_array2(classes_Train1,lon,lat,LON,LAT);
                fco2=field;
                biome_regimes1=classes_field;
                else
                fco2(1:180,1:360)=NaN;
                biome_regimes1(1:180,1:360)=NaN;
                end
            temp_array(month2plot,:,:)=fco2(:,:);
            biome_regimes2(month2plot,:,:)=biome_regimes1(:,:);
            end
                if(year2go==min(year_output))
                temp_array1=temp_array;
                biome_regimes3=biome_regimes2;
                else
                temp_array1=[temp_array1; temp_array];
                biome_regimes3=[biome_regimes3;biome_regimes2];
                end
        end

    %----------------------------------------------------------------------
    % if no observations per biome set to NAN
    %----------------------------------------------------------------------
    
    elseif(length(ind_FFN)<=1)
        hits_empty=hits_empty+1;
        hitlist(biome1)=3;
        temp_array1(1:max(year_output)-min(year_output)+1,1:180,1:360)=NaN;
        biome_regimes3(1:max(year_output)-min(year_output)+1,1:180,1:360)=NaN;
    end

%--------------------------------------------------------------------------
% 6) put the puzzle pieces together
%--------------------------------------------------------------------------

if(biome==content(1))
%if(biome==12)
    final_array=temp_array1;
    final_biomes=biome_regimes3;
else
    
    inddxx=find(isnan(final_array));
    final_array(inddxx)=temp_array1(inddxx);
    
    inddxx1=find(isnan(final_biomes));
    final_biomes(inddxx1)=biome_regimes3(inddxx1);
end

end

%fclose(fileID)8;
%--------------------------------------------------------------------------
% 7) create output array and save
%--------------------------------------------------------------------------

data_all=zeros(492,180,360);
data_all(:,:,:)=NaN;
data_biomes=zeros(492,180,360);
data_biomes(:,:,:)=NaN;
ind2=find(timevec==round(max(yearyeah+1)));
ind=find(timevec==round(min(yearyeah)));
data_all(ind:ind2-1,:,:)=final_array;
data_biomes(ind:ind2-1,:,:)=final_biomes;

clear ind ind2

saveloc=['output/BIOMEoutput_SOCAT/pCO2_Chl' net2take '_biome_' num2str(nnnumber) '.mat'];

save(saveloc, 'data_all', 'data_biomes','hitlist','hitdata');

hits_FFN
hits_mean
hits_empty

%end

%--------------------------------------------------------------------------
% 8) test plot of results
%--------------------------------------------------------------------------

xx=size(lat);
yy=size(lon);

latt=zeros(xx(2),yy(2));
lonn=zeros(xx(2),yy(2));
for ii=1:1:yy(2)
   latt(:,ii) = lat(:);   
end
for ii=1:1:xx(2)
   lonn(ii,:) = lon(:);  
end

plot_array1=squeeze(mode(data_biomes));
plot_array2=squeeze(nanmean(data_biomes));

for jj=1:180
    for kk=1:360
        ind=find(~isnan(data_biomes(:,jj,kk)));
        plot_array3(jj,kk)=numel(unique(data_biomes(ind,jj,kk)));
    end
end

plot_array3(plot_array3==0)=NaN;

map_plot_new(lon, lat, squeeze(nanmean(data_all(:,:,:))), 'equidistant', lonmin, lonmax, latmin, latmax, 1, 280, 450, 'fCO2',1,2)
map_plot_new(lon, lat, plot_array2, 'equidistant', lonmin, lonmax, latmin, latmax, 2, 1, biome, 'BGC mean',1,2)
%discrete_map_plot_long(lonn, latt, plot_array1, 'equidistant', lonmin, lonmax, latmin, latmax, 3, 1, biome, 'BGC mode')
%discrete_map_plot_long(lon, lat, plot_array3, 'equidistant', lonmin, lonmax, latmin, latmax, 4, 1, max(max(plot_array2)), ['number of regions per pixel with Longhurst'])

% time=min(year_output):1/12:max(year_output)+11/12;
% 
% hovmoeller_plot1(final_array,latt,time, 3,300, 420,'name')
% hovmoeller_plot1(final_biomes,latt,time, 4,1, biome,'name')