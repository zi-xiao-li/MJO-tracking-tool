

%MJO_tracking_2


%-----------------------------------------------------------------------------------

speed_slow = 3.3;%慢速MJO
speed_fast = 4.5;%快速MJO

file_name1 = "D:\project\output\PropagationInfo_fast.xlsx";
file_name2 = "D:\project\output\PropagationInfo_slow.xlsx";
file_name3 = "D:\project\output\PropagationInfo_all.xlsx";

%-----------------------------------------------------------------------------------

%读取数据
data_input  = [];
data        = [];
data_input  = xlsread('D:\project\output\PropagationInfo.xlsx');
data        = data_input;
%去除数据里的某些值
%aa1         = [];
%aa1         = find( data(:,10)>9 );
%data(aa1,:) = [];
aa1         = [];
aa1         = find(data(:,11)<=80 & data(:,12)>=120);
data        = data(aa1,:);
% slow fast 分类, 4m/s 和5m/s
speed     = [];
data_slow = [];
data_fast = [];
speed     = data(:,10);
data_slow = data(find(speed<speed_slow),:);  % 要特别注意这里是＜而不是≤
data_fast = data(find(speed>speed_fast),:);
nn_slow   = size(data_slow,1);
nn_fast   = size(data_fast,1);
xlswrite(file_name3,data)
xlswrite(file_name1,data_fast)
xlswrite(file_name2,data_slow)
% 计算对应速度
disp(['average speed of slow MJO:',num2str(mean(data_slow(:,10))),'m/s']);
disp(['average speed of fast MJO:',num2str(mean(data_fast(:,10))),'m/s']);
disp(['average speed of all MJO:',num2str(mean(data(:,10))),'m/s']);
disp(['fast MJO frequency:',num2str(size(data_fast,1))]);
disp(['slow MJO frequency:',num2str(size(data_slow,1))]);


    h = [];
    h = histogram(speed);
    h.FaceColor = [255 129 113]./255; % 颜色
    h.BinWidth = 0.5; % 设置直方图的宽度
    %h.EdgeColor = [255 129 113]./255;
    %坐标设置
    c              = colorbar;
    ax             = gca;
    title("a.                                                                               ",'FontSize',16)
    ylabel('Number of Events')
    xlabel('MJO Phase Speed (m/s)')
    c.LineWidth    = 1.;
    c.FontName     = 'Arial';
    c.FontSize     = 12;
    ax.YTick       = [0 2 4 6 8 10 12 14 16 18 20];
    ax.YTickLabel  = {'0', '2', '4', '6', '8', '10', '12', '14', '16', '18', ''};
    ax.YLim = [0 19];
    ax.FontSize    = 16;
    ax.LineWidth   = 2;
    ax.FontName    = 'Arial';
    position       = [0.140104166666667 0.179906542056075 0.362499999999999 0.74275700934577];
    ax.Position    = position;




%===================================================================================
%
%                                    合成分析
%
%===================================================================================

%all composite
All_MJO = [];
for ii=1:size(data,1)
    aa=[];
    aa=find(time(:,1)==data(ii,7) & time(:,2)==data(ii,8) &time(:,3)==data(ii,9));
    bb=[];
    bb=olr_EQ_ave(:,[aa-30:aa+30]);
    All_MJO(:,:,ii) = bb;
end
All_MJO_Composite = [];
for ii=1:size(All_MJO,1)
    for jj=1:size(All_MJO,2)
        All_MJO_Composite(ii,jj) = mean(All_MJO(ii,jj,:));
    end
end
%fast composite
Fast_MJO = [];
for ii=1:size(data_fast,1)
    aa=[];
    aa=find(time(:,1)==data_fast(ii,7) & time(:,2)==data_fast(ii,8) &time(:,3)==data_fast(ii,9));
    bb=[];
    bb=olr_EQ_ave(:,[aa-30:aa+30]);
    Fast_MJO(:,:,ii) = bb;
end
Fast_MJO_Composite = [];
for ii=1:size(Fast_MJO,1)
    for jj=1:size(Fast_MJO,2)
        Fast_MJO_Composite(ii,jj) = mean(Fast_MJO(ii,jj,:));
    end
end
%slow composite
Slow_MJO = [];
for ii=1:size(data_slow,1)
    aa=[];
    aa=find(time(:,1)==data_slow(ii,7) & time(:,2)==data_slow(ii,8) &time(:,3)==data_slow(ii,9));
    bb=[];
    bb=olr_EQ_ave(:,[aa-30:aa+30]);
    Slow_MJO(:,:,ii) = bb;
end
Slow_MJO_Composite = [];
for ii=1:size(Slow_MJO,1)
    for jj=1:size(Slow_MJO,2)
        Slow_MJO_Composite(ii,jj) = mean(Slow_MJO(ii,jj,:));
    end
end







%====================================================================
%           t test fast：其实这个合成应该是不需要t检验的
% 我觉得这个说法应该是不对的，因为chen都进行了检验
%====================================================================


%总体均值
Xn_slow  = [];
Xn_fast  = [];
std_slow = [];
std_fast = [];
t_slow   = [];
t_fast   = [];
n_slow   = size(Slow_MJO,3);
n_fast   = size(Fast_MJO,3);

for ii=1:size(All_MJO,1)
    for jj=1:size(All_MJO,2)
        Xn_slow(ii,jj)   = mean(Slow_MJO(ii,jj,:));
        Xn_fast(ii,jj)   = mean(Fast_MJO(ii,jj,:));
        std_slow(ii,jj)  = std(Slow_MJO(ii,jj,:));
        std_fast(ii,jj)  = std(Fast_MJO(ii,jj,:));
    end
end
% t
for ii=1:size(All_MJO,1)
    for jj=1:size(All_MJO,2)
        t_slow(ii,jj) = (Xn_slow(ii,jj))./(std_slow(ii,jj)./sqrt(n_slow));
        t_fast(ii,jj) = (Xn_fast(ii,jj))./(std_fast(ii,jj)./sqrt(n_fast));
    end
end









%====================================================================
%                               PLOT
%   修改一下，改成过显著性检验的contourf，然后画等值线
%====================================================================

%close all`

data_save = [];
t_save    = [];
data_save(:,:,1) = Fast_MJO_Composite([I6:I7],:)';
data_save(:,:,2) = Slow_MJO_Composite([I6:I7],:)';
position         = [0.140104166666667 0.179906542056075 0.362499999999999 0.74275700934577 ; 0.570340909090909 0.18107476635514 0.354659090909091 0.743925233644845];
t_save(:,:,1)    = t_fast([I6:I7],:)';
t_save(:,:,2)    = t_slow([I6:I7],:)';

% 把高于0.01显著性检验的提出来
% 尼玛，没用！
% aa1 = [];
% aa2 = [];
% aa1 = find(abs(t_save(:,:,1)) < 2.04); %n=31
% aa2 = find(abs(t_save(:,:,2)) < 2.05);
% tmp = data_save(:,:,1);     % 提取出该层
% tmp(aa1) = NaN;             % 在平面上按线性索引赋值
% data_save(:,:,1) = tmp;     % 再写回原数组
% tmp2 = data_save(:,:,2);
% tmp2(aa2) = NaN;
% data_save(:,:,2) = tmp2;


figure

%titile_str       = ["a. Fast              6.9m/s" "b. Slow              3.7m/s" ];
titile_str       = ["b. Fast (24 cases)                                           5.4 m/s" "c. Slow (25 cases)                                        2.7 m/s" ];

colortable = [0 23 69; 0 86 159; 0 113 184; 23 161 207; 145 220 231; 251 251 251;
251 251 251; 255 188 170; 255 106 70; 237 16 24; 194 6 12; 111 0 0]; 
colortable = colortable./255;
colormap(colortable);

%subplot(1,2,1)
%contourf(data_save(:,:,1));
%subplot(1,2,2)
%contourf(data_save(:,:,2));


for ii=1:2

    subplot(1,2,ii)
    %正值实线负值虚线
    z  = [];
    z1 = [];
    z2 = [];
    z  = data_save(:,:,ii);
    z1 = z;
    z2 = z;
    z1(z1>0)  = 0;
    z2(z2<=0) = 0;
    %-----------------------------------------------------------------------------------------
    contourf(z,'LineStyle','None');
    hold on
    contour(z1,'LineStyle','--','LineWidth',1.5,'LevelList',[-30:5:-5],'LineColor','k')
    hold on
    contour(z2,'LineStyle','-','LineWidth',1.5,'LevelList',[5:5:30],'LineColor','k')
    hold on
    %打点
    [XX YY] = meshgrid([1:size(data_save,2)],[1:size(data_save,1)]);
    FLAG = [];
    FLAG = find(abs(t_save(:,:,ii)) >= 2.75 );
    X=XX(FLAG);
    Y=YY(FLAG);
    %scatter(X,Y,5,'o','markerfacecolor','k','markeredgecolor','k')

    %-----------------------------------------------------------------------------------------
    %坐标设置
    c              = colorbar;
    ax             = gca;
    caxis([-30,30])
    set(c,'YTick',-30:5:30);
    title(titile_str(ii),'FontSize',16)
    ylabel('Days')
    c.LineWidth    = 1.;
    c.FontName     = 'Arial';
    c.FontSize     = 12;
    ax.FontSize     = 16;
    ax.LineWidth   = 2;
    ax.Position    = position(ii,:);
    X_Tick         = [5:12:81];
    ax.XTick       = X_Tick;
    X_Tick_Label_Pre = [];
    for i=I6:I7
        X_Tick_Label_Pre{i-I6+1} = Lon_Label{i};
    end
    X_Tick_Label =[];
    for i=1:length(X_Tick)
        X_Tick_Label{i}=X_Tick_Label_Pre{X_Tick(i)};
    end
    for i=1:length(X_Tick_Label)
        ax.XTickLabel(i)  = X_Tick_Label(i);
    end
    ax.YTick       = [1 11 21 31 41 51 61];
    ax.YTickLabel  = {'-30', '-20', '-10', '0', '10', '20', '30'};
    ax.FontName    = 'Arial';

    %-----------------------------------------------------------------------------------------
    %设置网格线
    %grid on
    %set(gca,'GridLineStyle','--','GridColor','k','GridAlpha',.4);  %GridAlpha 透明度
    %ax.GridLineWid = .4;

    

end

