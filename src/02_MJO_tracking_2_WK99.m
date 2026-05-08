%-------------------------------------------------------------------------
%
%  用于MJO tracking 截面的截取



clear all; close all; clc




%============================================================================
n        = 365;                        
t_start  = 1979;                       
t_end    = 2014;
lat_N    = 5;  % 5--5平均
lat_S    = -5;
IO_lon1  = 90;
IO_lon2  = 90; % 90°E为参考点
file_name = "D:\project\data\1\olr_MJO.nc";
%============================================================================

%给经度贴标签
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
% 计算符合clmDayTLL的日期，删掉2月29日
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
% 取赤道10--10平均
[~,I1]	= min(abs(lat(:)-lat_N));
[~,I2]	= min(abs(lat(:)-lat_S));
olr_EQ  = olr(:,[I1:I2],:);
for ii=1:size(olr_EQ,1)
    for jj=1:size(olr_EQ,3)
        olr_EQ_ave(ii,jj) = mean(olr_EQ(ii,:,jj));
    end
end
% 计算每一个longitude上的标准差SD
for ii=1:size(olr_EQ_ave,1)
    olr_std(ii,:)  = std(olr_EQ_ave(ii,:));
    olr_mean(ii,:) = mean(olr_EQ_ave(ii,:));
end
% 对MJO事件进行识别
% MJO index，90°E为参考点
[~,I4]       = min(abs(lon(:)-IO_lon1));
[~,I5]       = min(abs(lon(:)-IO_lon2));
for ii=1:size(olr_EQ,3)
    MJO_index(ii)  = mean(mean(olr_EQ([I4:I5],:,ii)));
end
olr_IO_std   = std(MJO_index);
olr_IO_mean  = mean(MJO_index);

% % 连续五天MJO index小于1SD,并选出五天中的最小值
% %! 说实话不知道你当初怎么想的，怎么会用这个方法？？
% MJO_date      = [];
% MJO_min_OLR   = [];
% for s = 1:length(MJO_index)-4
%     B = MJO_index(s:s+4); %获得连续5个数
%     if (length(find(B <= olr_IO_mean-olr_IO_std))==5) %判断5个元素是否大于等于IO_std
%         MJO_date    = [MJO_date s];
%         MJO_min_OLR = [MJO_min_OLR min(B)]; %找到五个数中最小值
%     end
%     %因为选取的范围是s:s+4，但只记录了s，所以必须把s+4还原回去，才是一次MJO事件的完整日期
%     %! 我明白了，你小子是把一次MJO定义为小于1SD的连续5天
%     %! 但我觉得这个定义有点问题，没有谁用这个定义啊？
%     %! 当然也不能这么说，neena也是把KW定义为三天，但我觉得你这个“完整日期”就有问题吧？我们的目的是要找到极小值日期
%     if length(MJO_date)>=2     
%         ll_1=MJO_date(end);
%         ll_2=MJO_date(end-1);
%         if ll_1-ll_2~=1
%             MJO_date = [MJO_date(1:end-1) ll_2+1 ll_2+2 ll_2+3 ll_2+4 MJO_date(end)];            
%         end
%         ll_3=MJO_min_OLR(end);
%         ll_4=MJO_min_OLR(end-1);
%         if ll_1-ll_2~=1
%             MJO_min_OLR = [MJO_min_OLR(1:end-1) MJO_min_OLR(end-1) MJO_min_OLR(end-1)...
%                     MJO_min_OLR(end-1) MJO_min_OLR(end-1) MJO_min_OLR(end)];            
%         end
%     end
% end

% MJO index极小值且<1SD
MJO_date     = [];
MJO_min_OLR  = [];
for s = 2:length(MJO_index)-1
    B = MJO_index(s);
    if B < MJO_index(s-1) & B < MJO_index(s+1) & B < olr_IO_mean-olr_IO_std 
        MJO_date = [MJO_date s];
        MJO_min_OLR = [MJO_min_OLR B];
    end
end
disp(['1979-2013年共发生MJO事件：',num2str(length(MJO_date)),'次'])

% % 冬季MJO事件选取,count1记录了冬季所有初始信息：日期，MJO index等
% count1 = [];
% count2 = [];
% var1   = [];
% var2   = [];
% var1       = time(MJO_date,:);
% var1(:,6)  = MJO_date;
% var1(:,7)  = MJO_min_OLR; % 全年
% var2       = find(var1(:,2)==11|var1(:,2)==12|var1(:,2)==1|...
%                         var1(:,2)==2|var1(:,2)==3|var1(:,2)==4); % 冬季
% count1     = var1(var2,:); %* 储存冬季变量
% % 确定起始日期和结束日期以及持续时间
% %! 说实话做这个有什么意义？？？找到最小值不就行了吗？？？
% var1=[];
% var2=[];
% var3=[];
% var1        = find(diff(count1(:,6))~=1);
% count_end   = count1(var1,:);             % ？？？？？？要不要加4
% count_end   = [count_end; count1(end,:)];
% count_start = count1(var1+1,:);
% count_start = [count1(1,:);count_start];
% count_day   = count_end(:,6)-count_start(:,6);
% disp(['1979-2013年冬季共发生MJO事件',num2str(length(count_day)),'次，其中：'])
% disp(['（1）主对流在IO区最长停留时间：',num2str(max(count_day+5)),'天'])
% disp(['（2）主对流在IO区最短停留时间：',num2str(min(count_day+5)),'天'])
% % 确定每次MJO事件对应的minist OLR index
% for ii=1:length(count_end)
%     var4          = [];
%     var4(:,[1:3]) = time([count_start(ii,6):count_end(ii,6)],[1:3]);
%     var4(:,4)     = MJO_index(count_start(ii,6):count_end(ii,6));
%     var2          = [var2 min(var4(:,4))]; %每次MJO中的lowest OLR
%     var3          = var4(find(var4(:,4)==min(var4(:,4))),[1:3]); % *lowest OLR对应日期
%     count2(ii,[1:3]) = var3;
% end
% count2(:,4) = var2; % 每次MJO事件中的lowest OLR日期及数值
% %! 你绕这么大一圈是要干嘛？？？？最后不还是得到OLR最小值对应的日期吗？

count2 = time(MJO_date,:);
count2 = count2(:,[1:3]); % 只保留日期
count2(:,4) = MJO_min_OLR; 

% 切片，为下一步做准备
[~,I6]       = min(abs(lon(:)-20));  % 切片在20°E-140°W
[~,I7]       = min(abs(lon(:)-220));
Segment_OLR_ALL = olr_EQ_ave([I6:I7],:);
Segment_Lon     = lon([I6:I7]);
Segment_t0      = count2(:,[1:3]); %! 服了，你看这儿，是不是就是minima的日期？？？
                                   %! 这个是极其重要的变量
% 截取1979.10.1-2014.5.31的数据
var1 = []; var2 = []; var3 = []; var4 = []; var5 = [];
var1 = find(time(:,1)==1979&time(:,2)==10&time(:,3)==1);  %!!!!!!!!!!!!!!!!!
var2 = find(time(:,1)==2014&time(:,2)==5&time(:,3)==31);  %!!!!!!!!!!!!!!!!!
var3 = Segment_OLR_ALL(:,[var1:var2]);
var4 = time([var1:var2],:);
% 截取boreal winter，注意这不是最终的boreal winter！
%? 是因为有5月吗？？
var5 = find(var4(:,2)==10|var4(:,2)==11|var4(:,2)==12|var4(:,2)==1|var4(:,2)==2|var4(:,2)==3|var4(:,2)==4|var4(:,2)==5);   %!!!!!!!!!!!!!!!!!
Segment_Time_Winter = var4(var5,:);
Segment_OLR_Winter  = var3(:,var5);
% 循环得到每一年冬季的hov图切片segment
Segment_anual_number    = 243; %!!!!!!!!!!!!!!!!! 每一年冬季（包含5月）的天数
Segment_year_number     = size(Segment_Time_Winter,1)/Segment_anual_number;
for ii=1:Segment_year_number
    Segment_OLR(:,:,ii) = Segment_OLR_Winter(:,[1+(ii-1)*Segment_anual_number:...
                                    ii*Segment_anual_number]);
    Segment_Time(:,:,ii)= Segment_Time_Winter([1+(ii-1)*Segment_anual_number:...
                                    ii*Segment_anual_number],:);
end

%把每一年冬季的参考日期提出来
%*确实是这样的，就是把每一年冬季的Segment_t0提出来
%!! 这里不太对吧？？？？上面切片用5月10月还可以理解，这里没必要吧？？？
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



