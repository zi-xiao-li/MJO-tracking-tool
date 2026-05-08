
%MJO_tracking_2

%————————————————————————————————————————————————————————————————————————————

data_print=[];

for Z=1:35

The_Chosen_One = 1979+Z-1;

%————————————————————————————————————————————————————————————————————————————

ex_time    = []; 
ex_OLR     = []; 
ex_t0      = []; 
ex_lon     = []; 
ex_y0      = []; 
ex_x0      = [];
EX_t0      = [];
ex_ref_y0  = [];

Propagate_Start     = [];
Propagate_End       = [];
Propagate_Day_0     = [];
Propagate_Start_Lon = [];
Propagate_End_Lon   = [];
Propagate_Speed     = [];
Propagate_Bm        = [];

%找到年份对应的序号
Number_All = [1979:2013];
Number     = find(Number_All==The_Chosen_One);

ex_OLR     = Segment_OLR(:,:,Number)';
ex_x0      = 90;   % 对应经度
ex_time    = Segment_Time(:,[1:3],Number);
ex_lon     = Segment_Lon;
EX_t0      = Segment_t0{Number}; %即2006年的t0,以年月日的形式

%如果该年没有发生MJO事件>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
if isempty(EX_t0)==1
    continue
end

















for H=1:size(EX_t0,1)



    ex_t0(1,:) = EX_t0(H,:);

% 在181day的hov图中y轴：t∈[t0-12:t0+12],t0的坐标
ex_ref_y0  = find(ex_time(:,1)==ex_t0(1,1) & ex_time(:,2)==ex_t0(1,2) & ex_time(:,3)==ex_t0(1,3));
% % 跳过越界事件，比如5.29
% % *这个其实是合理的，因为我们选择的日期是10，11，12，1，2，3，4，5月，太远的应该没关系
% % *大不了之后再删了就行
%!之前出现越界问题是因为Segment_t0选取了5,10月，删除后不会出现这个问题
%!本身切片选择5，10月就是为了提前解决越界可能性
% if ex_ref_y0 - 12 < 1 | ex_ref_y0 + 12 > size(ex_OLR, 1)
%     disp(['跳过越界事件：', num2str(ex_t0(1,:))])
%     continue
% end
ex_y0      = [ex_ref_y0-12:ex_ref_y0+12]-0.5; % ref的y0前后12d



%如果t∈[t0-12,t0+12]中有小于0的数，比如参考点在11.1日>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%if isempty(find(ex_y0<0))==0
%    ex_y0 = ex_y0(find(ex_y0>0)); %因为后面要yo-0.5,x0-0.5，所以0的那项去掉。得到y轴坐标
%end



%如果t∈[t0-12,t0+12]中有大于181的数，比如参考点在4.30日>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%if isempty(find(ex_y0>181))==0
%    ex_y0 = ex_y0(find(ex_y0<181)); %因为后面要yo-0.5,x0-0.5，所以0的那项去掉。得到y轴坐标
%end




% 通过t∈[y0-12,y0+12]，找到直线
syms x y

[~,x0]    = min(abs(ex_lon(:)-ex_x0));%参考点横坐标是不变的，就是90°E
x0        = x0 - 0.5;
%计算K值，每隔0.1m/s
K_real = [1:0.1:25];
K      = 111001./(34560.*K_real+1);


%从这里开始在坐标轴中表示>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

% 确定右半部分每一条直线与网格的交点
coords_left  = [];
coords_right = [];

coords_right = zeros(length(ex_y0),length(K),2); % 区域日期内过参考点的每个斜率与网格边缘的交点
for ii=1:length(ex_y0)
    y0 = ex_y0(ii) ;
    for jj=1:length(K)
        k  = K(jj);
        y = @(x) k*(x-x0) + y0;
        if (y(81) <= 243) % 如果斜率没有超过网格，则坐标为(81,y(81))
            coords_right(ii,jj,1) = 81;
            coords_right(ii,jj,2) = y(81);
        else 
            y = @(x) (1./k)*(x-y0)+x0; % 如果超过了，就为(x(181),181)
            coords_right(ii,jj,1) = y(243);
            coords_right(ii,jj,2) = 243;
        end
    end
end

% 确定左半部分每一条直线与网格的交点
coords_left = zeros(length(ex_y0),length(K),2); 
for ii=1:length(ex_y0)
    y0 = ex_y0(ii) ;
    for jj=1:length(K)
        k  = K(jj);
        y = @(x) k*(x-x0) + y0;
        if (y(0) >= 0) % 如果斜率没有超过网格，则坐标为(0,y(0))
            coords_left(ii,jj,1) = 0;
            coords_left(ii,jj,2) = y(0);
        else 
            y = @(x) (1./k)*(x-y0)+x0; % 如果超过了，就为(x(0),0)
            coords_left(ii,jj,1) = y(0);
            coords_left(ii,jj,2) = 0;
        end
    end
end

% 求坐标点连线经过的区域，即xxx日，斜率为xxx时已知左右顶点，求连线经过的区域
segs_info = [];
segs_info = cell(size(coords_left,1),size(coords_left,2));
for ii=1:size(coords_right,1)
    y0 = ex_y0(ii);
    for jj=1:size(coords_right,2)
        k     =[];
        P1    =[];
        P2    =[];
        xdx   =[];
        xdy   =[];
        ydx   =[];
        ydy   =[];
        k        = K(jj);
        % 找到经确认的起始点
        P1(:,:)  = coords_left(ii,jj,:);
        P2(:,:)  = coords_right(ii,jj,:);
        P1_x     = P1(1);
        P1_y     = P1(2);
        P2_x     = P2(1);
        P2_y     = P2(2);
        xmin     = min(P1_x,P2_x);
        ymin     = min(P1_y,P2_y);
        xmax     = max(P1_x,P2_x);
        ymax     = max(P1_y,P2_y);
        % x取整对应的y坐标    ???????????????????????????????????????????
        y = @(x) k*(x-x0) + y0;
        for i = ceil(xmin):floor(xmax)
            xdx(i-ceil(xmin)+1) = i;
            xdy(i-ceil(xmin)+1) = y(i);
            if (y(i)<0 && y(i)>-1)
                xdy(i-ceil(xmin)+1) = 0;
                disp(['有数字小于0，具体数值为',num2str(y(i)),'参考点：',num2str(x0),',',num2str(y0),...
                        '位置：',num2str(ii),'天，','斜率为',num2str(jj),'xy轴坐标为',num2str(i),',',num2str(y(i))]);
            elseif (y(i)<-1)
                disp(['有数字小于-1，具体数值为',num2str(y(i)),'参考点：',num2str(x0),',',num2str(y0),...
                        '位置：',num2str(ii),'天，','斜率为',num2str(jj),'xy轴坐标为',num2str(i),',',num2str(y(i))]);
            end
        end
        % y取整对应的x坐标
        y = @(x) (1./k)*(x-y0)+x0;
        for j = ceil(ymin):floor(ymax)
            ydy(j-ceil(ymin)+1) = j;
            ydx(j-ceil(ymin)+1) = y(j);
            if (y(j)<0 && y(j)>-1)
                ydx(j-ceil(ymin)+1) = 0;
                disp(['有数字小于0，具体数值为',num2str(y(j)),'参考点：',num2str(x0),',',num2str(y0),...
                        '位置：',num2str(ii),'天，','斜率为',num2str(jj),'xy轴坐标为',num2str(y(j)),',',num2str(j)]);
            elseif (y(j)<-1)
                disp(['有数字小于-1，具体数值为',num2str(y(j)),'参考点：',num2str(x0),',',num2str(y0),...
                        '位置：',num2str(ii),'天，','斜率为',num2str(jj),'xy轴坐标为',num2str(y(j)),',',num2str(j)]);
            end
        end

        SP = [];
        SP = unique([P1_x,P2_x,xdx,ydx; P1_y,P2_y,xdy,ydy]','rows'); %单调增函数确实可以这么做

        segs = [];
        for t = 1 : size(SP,1)-1
            segs(t).index_x = max(ceil(SP(t+1,1)),ceil(SP(t,1))); % 列
            segs(t).index_y = max(ceil(SP(t+1,2)),ceil(SP(t,2))); % 行
            segs(t).OLRA    = ex_OLR(segs(t).index_y, segs(t).index_x);
        end
        segs_info{ii,jj} = segs; %得到的是直线穿过的格点序号
    end
end

% 找到每条trail line OLRA小于对应平均值以下一个标准差的最长片段，longitude gap < 10°仍为一个片段;
% 我们只在意每个整数经度上的值
for ii=1:size(segs_info,1)
    for jj=1:size(segs_info,2)
        var1   = [];
        Series = [];
        var1   = segs_info{ii,jj};

        %识别所有<平均值下1SD的格点
        for t=1:length(var1)
            if var1(t).OLRA <= Segment_mean(find(Segment_std(:,2)==var1(t).index_x)) ...
                                -Segment_std(find(Segment_std(:,2)==var1(t).index_x))
            %if var1(t).OLRA <= -10
                var1(t).track = 1;
            else 
                var1(t).track = 0;
            end
            Series(t)=var1(t).track;
        end

        % 找到0数据段，如果相隔不超过10°，则仍为一次事件
        temp = [];aa_start=[]; aa_end=[];
        temp = diff(Series);
        aa_start=find(temp==-1)+1; %0序列开始
        aa_end=find(temp==1);  %0序列的结束

        if Series(1)==0
            aa_start=[1 aa_start];
        end

        if Series(end)==0
            aa_end=[aa_end length(Series)];
        end

        for i=1:length(aa_end)
            if var1(aa_end(i)).index_x - var1(aa_start(i)).index_x <=4
                Series(aa_start(i):aa_end(i)) = 1;
            end
        end

        for t=1:length(segs_info{ii,jj})
            segs_info{ii,jj}(t).track = Series(t);
        end

        % Longest Segment，找到1序列中最长的片段，最长指的是不要那些细碎片段
        temp = []; aa_start=[]; aa_end=[];
        temp = diff(Series);
        aa_start=find(temp==1)+1; %1序列开始
        aa_end=find(temp==-1);  %1序列的结束

        if Series(1)==1
            aa_start=[1 aa_start];
        end

        if Series(end)==1
            aa_end=[aa_end length(Series)];
        end

        % 计算A(t,c),L(t,c)
        a=[]; b=[]; bb_start=[]; bb_end=[]; cc=0;
        if isempty(aa_end)==1
            segs_info{ii,jj}(1).StartLong=0;
            segs_info{ii,jj}(1).EndtLong=0;
            segs_info{ii,jj}(1).A = 0;
            segs_info{ii,jj}(1).L = 0;
            segs_info{ii,jj}(1).Start=0;
            segs_info{ii,jj}(1).End=0;
        else
            [a b]=max(aa_end-aa_start);
            bb_start = aa_start(b);
            bb_end   = aa_end(b); 
            for t=bb_start:bb_end
                 cc=cc+segs_info{ii,jj}(t).OLRA;
            end
            segs_info{ii,jj}(1).StartLong=segs_info{ii,jj}(bb_start).index_x;
            segs_info{ii,jj}(1).EndLong=segs_info{ii,jj}(bb_end).index_x;
            segs_info{ii,jj}(1).StartDate=segs_info{ii,jj}(bb_start).index_y;
            segs_info{ii,jj}(1).EndDate=segs_info{ii,jj}(bb_end).index_y;
            segs_info{ii,jj}(1).A=cc;
            segs_info{ii,jj}(1).L=2.5*(segs_info{ii,jj}(bb_end).index_x - segs_info{ii,jj}(bb_start).index_x);
        end
    end
end

A = [];
L = [];
B = [];
A_m = [];
L_m = [];
B_m = [];
pos_x = [];
pos_y = [];
for ii=1:size(segs_info,1)
    for jj=1:size(segs_info,2)
        A(ii,jj) = abs(segs_info{ii,jj}(1).A);
        L(ii,jj) = abs(segs_info{ii,jj}(1).L);
    end
end

%找到Am、Lm
A_m = max(max(A));
L_m =max(max(L));
%B(t,c)
B = A./A_m + L./L_m;
B_m = max(max(B));
[pos_x pos_y] = find(B==B_m);
%对应的日期和斜率,posx日期，posy斜率

%K斜率间隔过小可能导致出现多个pos_y的问题，不仅在斜率，天数也一样>>>>>>>>>>>>>>>>>>>>>>>>>>>

if length(pos_x)~=1
    disp(['pos x(日期)',num2str(pos_x')])
    pos_x = max(pos_x); %可能导致生命周期比较长？
    %pos_x = nan;
end

if length(pos_y)~=1
    disp(['pos y(斜率)',num2str(pos_y')])
    pos_y = max(pos_y);  %速度偏慢
    %pos_y = nan;
end

Final_K    = [];
Final_info = [];
Final_K    = K_real(pos_y);
Final_info = segs_info{pos_x,pos_y};


disp(['Start：',num2str(ex_time(Final_info(1).StartDate,:))])
disp(['End：',num2str(ex_time(Final_info(1).EndDate,:))])
disp(['Day 0：',num2str(ex_t0)])
disp(['Start Lon：',num2str(Final_info(1).StartLong*2.5 + 20)])
disp(['End Lon：',num2str(Final_info(1).EndLong*2.5 + 20)])
disp(['MJO传播速度：',num2str(Final_K),'m/s'])

Propagate_Start     = [Propagate_Start; ex_time(Final_info(1).StartDate,:)];
Propagate_End       = [Propagate_End; ex_time(Final_info(1).EndDate,:)];
Propagate_Day_0     = [Propagate_Day_0; ex_t0];
Propagate_Start_Lon = [Propagate_Start_Lon; Final_info(1).StartLong*2.5 + 20];
Propagate_End_Lon   = [Propagate_End_Lon; Final_info(1).EndLong*2.5 + 20];
Propagate_Speed     = [Propagate_Speed; Final_K];
Propagate_Bm        = [Propagate_Bm; B_m];

end

data=[];
data=[Propagate_Start Propagate_End Propagate_Day_0 Propagate_Speed Propagate_Start_Lon Propagate_End_Lon Propagate_Bm];
data_print=[data_print;data];


end




file_name1 = "D:\project\output\PropagationInfo.xlsx";  %% 注意把重复值剔除！以及日期不对的
file_name2 = "D:\project\output\PropagationInfo_prepare.xlsx";

xlswrite(file_name1,data_print)

%读取数据
data        = [];
data        = xlsread(file_name1);
%去除数据里的某些值
aa1         = [];
aa1         = find(data(:,11)<=80 & data(:,12)>=120);
data        = data(aa1,:);

xlswrite(file_name2,data)

