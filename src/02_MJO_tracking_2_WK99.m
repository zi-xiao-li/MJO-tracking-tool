%-----------------------------------------------------------------------
% Used for extracting sections for MJO tracking
%-----------------------------------------------------------------------


clear all; close all; clc




%============================================================================
n        = 365;                        
t_start  = 1979;                       
t_end    = 2014;
lat_N    = 5;  % 5x5 averaging  
lat_S    = -5;
IO_lon1  = 90;
IO_lon2  = 90; % 90°E as the reference point 
file_name = "D:\project\data\1\olr_MJO.nc";
%============================================================================

% Label the longitudes 
aa        = [];
aa        = [0:2.5:180];
Lon_Label = [];
Lon_Label{1} = num2str(0);
for ii=2:length(aa)-1
    Lon_Label{ii} = strcat(num2str((aa(ii))),'°E');
end
Lon_Label{length(aa)} = num2str(180);
aa =[];
aa = [182.5:2.5:357.5];
for ii=1:length(aa)
    Lon_Label{73+ii} = strcat(num2str(360-aa(ii)),'°W');
end

info     = ncinfo(file_name);
lat      = ncread(file_name,'lat');
lon      = ncread(file_name,'lon');
olr      = ncread(file_name,'olr');
time_pre = ncread('D:\project\data\0\olr.day.mean.nc','time');
time_pre = (time_pre-time_pre(1))./24.;

% transfer to Julian day
t1       = datetime('1974-06-01');
jd1      = juliandate(t1);
time_pre = time_pre+jd1;
for ii=1:length(time_pre);
    [a b c]    = julian2greg(time_pre(ii));
	time(ii,1) = a;
	time(ii,2) = b;
	time(ii,3) = c;
end
time_pre = time;
a=[];
b=[];
c=[];
% cut off the time,1979-2013
a     = find(time(:,1)==t_start&time(:,2)==1&time(:,3)==1);
b     = find(time(:,1)==t_end&time(:,2)==12&time(:,3)==31);
time  = time_pre([a:b],:,:);
% Calculate dates matching clmDayTLL, removing February 29 
for ii=1:t_end-t_start+1
	c         = find(time(:,1)==t_start+ii-1);
	d         = length(c);
	time(c,4) = [1:d];
end
time(:,5)     = time(:,1).*1000+time(:,4);
a             = find(time(:,2)==2&time(:,3)==29);
time(a,:)     = [];
for ii=1:length(find(time(:,4)==365))
    b(n*(ii-1)+1:n*ii) = linspace(1,n,n);
end
time(:,4) = b;
% Take 10°S–10°N zonal mean  
[~,I1]	= min(abs(lat(:)-lat_N));
[~,I2]	= min(abs(lat(:)-lat_S));
olr_EQ  = olr(:,[I1:I2],:);
for ii=1:size(olr_EQ,1)
    for jj=1:size(olr_EQ,3)
        olr_EQ_ave(ii,jj) = mean(olr_EQ(ii,:,jj));
    end
end
% Calculate standard deviation (SD) at each longitude  
for ii=1:size(olr_EQ_ave,1)
    olr_std(ii,:)  = std(olr_EQ_ave(ii,:));
    olr_mean(ii,:) = mean(olr_EQ_ave(ii,:));
end
% Identify MJO events  
% MJO index, 90°E as reference point 
[~,I4]       = min(abs(lon(:)-IO_lon1));
[~,I5]       = min(abs(lon(:)-IO_lon2));
for ii=1:size(olr_EQ,3)
    MJO_index(ii)  = mean(mean(olr_EQ([I4:I5],:,ii)));
end
olr_IO_std   = std(MJO_index);
olr_IO_mean  = mean(MJO_index);

% MJO index minimum value and < 1 SD
MJO_date     = [];
MJO_min_OLR  = [];
for s = 2:length(MJO_index)-1
    B = MJO_index(s);
    if B < MJO_index(s-1) & B < MJO_index(s+1) & B < olr_IO_mean-olr_IO_std 
        MJO_date = [MJO_date s];
        MJO_min_OLR = [MJO_min_OLR B];
    end
end
disp(['Total MJO events from 1979 to 2013: ', num2str(length(MJO_date))])


count2 = time(MJO_date,:);
count2 = count2(:,[1:3]); 
count2(:,4) = MJO_min_OLR; 

% Slice data, preparing for the next step
[~,I6]       = min(abs(lon(:)-20));  % Slices cover 20°E–140°W
[~,I7]       = min(abs(lon(:)-220));
Segment_OLR_ALL = olr_EQ_ave([I6:I7],:);
Segment_Lon     = lon([I6:I7]);
Segment_t0      = count2(:,[1:3]); 
                                   
% Select data from 1979-10-01 to 2014-05-31
var1 = []; var2 = []; var3 = []; var4 = []; var5 = [];
var1 = find(time(:,1)==1979&time(:,2)==10&time(:,3)==1);  %!!!!!!!!!!!!!!!!!
var2 = find(time(:,1)==2014&time(:,2)==5&time(:,3)==31);  %!!!!!!!!!!!!!!!!!
var3 = Segment_OLR_ALL(:,[var1:var2]);
var4 = time([var1:var2],:);
% Select boreal winter (note: this is not the final boreal winter!)
var5 = find(var4(:,2)==10|var4(:,2)==11|var4(:,2)==12|var4(:,2)==1|var4(:,2)==2|var4(:,2)==3|var4(:,2)==4|var4(:,2)==5);   %!!!!!!!!!!!!!!!!!
Segment_Time_Winter = var4(var5,:);
Segment_OLR_Winter  = var3(:,var5);
% Loop to obtain Hovmöller plot segments for each winter
Segment_anual_number    = 243;% Total days in each winter (including May)
Segment_year_number     = size(Segment_Time_Winter,1)/Segment_anual_number;
for ii=1:Segment_year_number
    Segment_OLR(:,:,ii) = Segment_OLR_Winter(:,[1+(ii-1)*Segment_anual_number:...
                                    ii*Segment_anual_number]);
    Segment_Time(:,:,ii)= Segment_Time_Winter([1+(ii-1)*Segment_anual_number:...
                                    ii*Segment_anual_number],:);
end

% Extract reference dates for each winter
% Extract Segment_t0 for each winter
Segment_t0_seperate = {};
for ii=1:Segment_year_number
    aaa1=[];
    aaa2=[];
    aaa1=Segment_t0(find(Segment_t0(:,1)==1978+ii),:);
    aaa2=Segment_t0(find(Segment_t0(:,1)==1979+ii),:);
    % First_Year=aaa1(find(aaa1(:,2)==10|aaa1(:,2)==11| aaa1(:,2)==12),:);
    % Last_Year =aaa2(find(aaa2(:,2)==1| aaa2(:,2)==2| aaa2(:,2)==3| aaa2(:,2)==4| aaa2(:,2)==5),:);
    First_Year=aaa1(find(aaa1(:,2)==11| aaa1(:,2)==12),:);
    Last_Year =aaa2(find(aaa2(:,2)==1| aaa2(:,2)==2| aaa2(:,2)==3| aaa2(:,2)==4),:);
    Segment_t0_seperate{ii} = [First_Year;Last_Year];
end

Segment_t0 = [];
Segment_t0 = Segment_t0_seperate;
Segment_std      = olr_std([I6:I7]);
Segment_std(:,2) = Segment_Lon;
Segment_mean     = olr_mean([I6:I7]);

clear Segment_Time_Winter Segment_OLR_Winter Segment_t0_seperate



