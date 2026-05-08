
%-----------------------------------------------------------------------
% 1. box的选取，进行数据预处理，后续进行wkSpaceTime处理
% 2. 预处理：1979-2013年OLR数据，去掉2月29日（去掉闰年，365对齐）
% 3. 虽然2013冬季是到14年4月，但是选box的数据是到2013年12月31日，参考WK99
%-----------------------------------------------------------------------     


clear all; close all; clc


fn      = 'D:\project\data\0\olr.day.mean.nc';
n       = 365;                        
t_start = 1979;                       
t_end   = 2013;


info     = ncinfo(fn);
lat      = ncread(fn,'lat');
lon      = ncread(fn,'lon');
olr_pre  = ncread(fn,'olr');
time_pre = ncread(fn,'time');
time_pre = (time_pre-time_pre(1))./24.;

% transfer to Julian day
t1       = datetime('1974-06-01');  % Daily values for 1974/06 - 2022/12/31
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
% cut off the time
a     = find(time(:,1)==t_start&time(:,2)==1&time(:,3)==1);
b     = find(time(:,1)==t_end&time(:,2)==12&time(:,3)==31);
olr   = olr_pre(:,:,[a:b]);
time  = time_pre([a:b],:,:);
% 计算符合clmDayTLL的日期，删掉2月29日
% 主要目的是把闰年去掉，便于计算吧，365d一个维度
for ii=1:t_end-t_start+1
	c         = find(time(:,1)==t_start+ii-1);
	d         = length(c);
	time(c,4) = [1:d];
end
time(:,5)   = time(:,1).*1000+time(:,4);
a           = find(time(:,2)==2&time(:,3)==29);
time(a,:)   = [];
olr(:,:,a)  = [];
a           = length(find(time(:,4)==365));
for ii=1:a
    b(n*(ii-1)+1:n*ii)=linspace(1,n,n);
end
time(:,4) = b;



clear ncid

ncid=netcdf.create('D:\project\data\1\olr_1979_2013.nc','CLOBBER');
%% 2 定义维度
dimidx = netcdf.defDim(ncid,'lon',size(olr,1)); 
dimidy = netcdf.defDim(ncid,'lat',size(olr,2));    
dimidz = netcdf.defDim(ncid,'time',size(olr,3));
%% 3 赋予属性
%lon
varid_lon=netcdf.defVar(ncid,'lon','double',dimidx);
netcdf.putAtt(ncid,varid_lon,'standard_name','Longitude');
netcdf.putAtt(ncid,varid_lon,'long_name','Longitude');
netcdf.putAtt(ncid,varid_lon,'units','degrees_east');
netcdf.putAtt(ncid,varid_lon,'axis','X');
%lat
varid_lat=netcdf.defVar(ncid,'lat','double',dimidy);
netcdf.putAtt(ncid,varid_lat,'standard_name','Latitude');
netcdf.putAtt(ncid,varid_lat,'long_name','Latitude');
netcdf.putAtt(ncid,varid_lat,'units','degrees_north');
%time
varid_time=netcdf.defVar(ncid,'time','double',dimidz);
netcdf.putAtt(ncid,varid_time,'standard_name','Time');
netcdf.putAtt(ncid,varid_time,'calendar','standard');
netcdf.putAtt(ncid,varid_time,'units','days since 1979-01-01 00:00:00');
%OLRA
varid_olr=netcdf.defVar(ncid,'olr','double',[dimidx dimidy dimidz]);
netcdf.putAtt(ncid,varid_olr,'_FillValue',-9999);
netcdf.putAtt(ncid,varid_olr,'missing_value',-9999);
%% 4 完成netCDF文件定义模式
netcdf.endDef(ncid)
%% 5 把数据写到netcdf的文件中
netcdf.putVar(ncid,varid_olr,olr);
netcdf.putVar(ncid,varid_lon,lon);
netcdf.putVar(ncid,varid_lat,lat);
netcdf.putVar(ncid,varid_time,[1:size(olr,3)]);  
%% 6 关闭文件
netcdf.close(ncid);
info_box = ncinfo('D:\project\data\1\olr_1979_2013.nc');

% 接下来衔接spectral analysis.ncl
