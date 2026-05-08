% Plot Hovmöller diagrams

% MJO_tracking_2

%-----------------------------------------------------------------------------------
speed_slow = 3.3; % Slow MJO
speed_fast = 4.5; % Fast MJO

file_name1 = "D:\project\output\PropagationInfo_fast.xlsx";
file_name2 = "D:\project\output\PropagationInfo_slow.xlsx";
file_name3 = "D:\project\output\PropagationInfo_all.xlsx";

%-----------------------------------------------------------------------------------
% Read data
data_input = [];
data       = [];
data_input = xlsread('D:\project\output\PropagationInfo.xlsx');
data       = data_input;

% Remove certain data points based on longitude
aa1       = [];
aa1       = find(data(:,11) <= 80 & data(:,12) >= 120);
data      = data(aa1,:);

% Classify slow and fast MJO (thresholds 4 m/s and 5 m/s)
speed      = data(:,10);
data_slow  = data(speed < speed_slow, :);  % Note: <, not ≤
data_fast  = data(speed > speed_fast, :);
nn_slow    = size(data_slow,1);
nn_fast    = size(data_fast,1);

% Write data to Excel
xlswrite(file_name3, data)
xlswrite(file_name1, data_fast)
xlswrite(file_name2, data_slow)

% Compute average speeds
disp(['Average speed of slow MJO: ', num2str(mean(data_slow(:,10))), ' m/s']);
disp(['Average speed of fast MJO: ', num2str(mean(data_fast(:,10))), ' m/s']);
disp(['Average speed of all MJO: ', num2str(mean(data(:,10))), ' m/s']);
disp(['Fast MJO frequency: ', num2str(size(data_fast,1))]);
disp(['Slow MJO frequency: ', num2str(size(data_slow,1))]);

% Plot histogram
h = histogram(speed);
h.FaceColor = [255 129 113]./255; % Set color
h.BinWidth  = 0.5; % Bin width
c           = colorbar;
ax          = gca;
title("a.                                                                               ",'FontSize',16)
ylabel('Number of Events')
xlabel('MJO Phase Speed (m/s)')
c.LineWidth = 1.;
c.FontName  = 'Arial';
c.FontSize  = 12;
ax.YTick    = [0 2 4 6 8 10 12 14 16 18 20];
ax.YTickLabel = {'0','2','4','6','8','10','12','14','16','18',''};
ax.YLim     = [0 19];
ax.FontSize = 16;
ax.LineWidth = 2;
ax.FontName = 'Arial';
ax.Position = [0.140104166666667 0.179906542056075 0.362499999999999 0.74275700934577];

%===================================================================================
%                             Composite Analysis
%===================================================================================

% All composite
All_MJO = [];
for ii = 1:size(data,1)
    aa = find(time(:,1) == data(ii,7) & time(:,2) == data(ii,8) & time(:,3) == data(ii,9));
    bb = olr_EQ_ave(:, [aa-30:aa+30]);
    All_MJO(:,:,ii) = bb;
end

% Compute mean for all MJO composite
All_MJO_Composite = [];
for ii = 1:size(All_MJO,1)
    for jj = 1:size(All_MJO,2)
        All_MJO_Composite(ii,jj) = mean(All_MJO(ii,jj,:));
    end
end

% Fast composite
Fast_MJO = [];
for ii = 1:size(data_fast,1)
    aa = find(time(:,1) == data_fast(ii,7) & time(:,2) == data_fast(ii,8) & time(:,3) == data_fast(ii,9));
    bb = olr_EQ_ave(:, [aa-30:aa+30]);
    Fast_MJO(:,:,ii) = bb;
end

Fast_MJO_Composite = [];
for ii = 1:size(Fast_MJO,1)
    for jj = 1:size(Fast_MJO,2)
        Fast_MJO_Composite(ii,jj) = mean(Fast_MJO(ii,jj,:));
    end
end

% Slow composite
Slow_MJO = [];
for ii = 1:size(data_slow,1)
    aa = find(time(:,1) == data_slow(ii,7) & time(:,2) == data_slow(ii,8) & time(:,3) == data_slow(ii,9));
    bb = olr_EQ_ave(:, [aa-30:aa+30]);
    Slow_MJO(:,:,ii) = bb;
end

Slow_MJO_Composite = [];
for ii = 1:size(Slow_MJO,1)
    for jj = 1:size(Slow_MJO,2)
        Slow_MJO_Composite(ii,jj) = mean(Slow_MJO(ii,jj,:));
    end
end

%====================================================================
%             t-test for fast/slow composites
%====================================================================
Xn_slow  = [];
Xn_fast  = [];
std_slow = [];
std_fast = [];
t_slow   = [];
t_fast   = [];
n_slow   = size(Slow_MJO,3);
n_fast   = size(Fast_MJO,3);

for ii = 1:size(All_MJO,1)
    for jj = 1:size(All_MJO,2)
        Xn_slow(ii,jj)  = mean(Slow_MJO(ii,jj,:));
        Xn_fast(ii,jj)  = mean(Fast_MJO(ii,jj,:));
        std_slow(ii,jj) = std(Slow_MJO(ii,jj,:));
        std_fast(ii,jj) = std(Fast_MJO(ii,jj,:));
    end
end

% Compute t-values
for ii = 1:size(All_MJO,1)
    for jj = 1:size(All_MJO,2)
        t_slow(ii,jj) = Xn_slow(ii,jj) / (std_slow(ii,jj)/sqrt(n_slow));
        t_fast(ii,jj) = Xn_fast(ii,jj) / (std_fast(ii,jj)/sqrt(n_fast));
    end
end

%====================================================================
%                               Plotting
%   Use contourf with significance overlay, then plot contours
%====================================================================
data_save = [];
t_save    = [];
data_save(:,:,1) = Fast_MJO_Composite([I6:I7],:)';
data_save(:,:,2) = Slow_MJO_Composite([I6:I7],:)';
position         = [0.140104166666667 0.179906542056075 0.362499999999999 0.74275700934577;
                    0.570340909090909 0.18107476635514 0.354659090909091 0.743925233644845];
t_save(:,:,1)    = t_fast([I6:I7],:)';
t_save(:,:,2)    = t_slow([I6:I7],:)';

% Define color map
titile_str = ["b. Fast (24 cases)                                           5.4 m/s", ...
              "c. Slow (25 cases)                                        2.7 m/s"];
colortable = [0 23 69; 0 86 159; 0 113 184; 23 161 207; 145 220 231; 251 251 251;
              251 251 251; 255 188 170; 255 106 70; 237 16 24; 194 6 12; 111 0 0] ./ 255;
colormap(colortable);

for ii = 1:2
    subplot(1,2,ii)
    z  = data_save(:,:,ii);
    z1 = z; z2 = z;
    z1(z1>0)  = 0;
    z2(z2<=0) = 0;

    % Filled contour for data
    contourf(z,'LineStyle','None'); hold on;
    % Contour for negative values
    contour(z1,'LineStyle','--','LineWidth',1.5,'LevelList',[-30:5:-5],'LineColor','k'); hold on;
    % Contour for positive values
    contour(z2,'LineStyle','-','LineWidth',1.5,'LevelList',[5:5:30],'LineColor','k'); hold on;

    % Overlay significant points
    [XX, YY] = meshgrid(1:size(data_save,2), 1:size(data_save,1));
    FLAG = find(abs(t_save(:,:,ii)) >= 2.75);
    X = XX(FLAG);
    Y = YY(FLAG);
    %scatter(X,Y,5,'o','markerfacecolor','k','markeredgecolor','k')

    % Axis settings
    c           = colorbar;
    caxis([-30,30])
    set(c,'YTick',-30:5:30);
    title(titile_str(ii),'FontSize',16)
    ylabel('Days')
    c.LineWidth = 1.;
    c.FontName  = 'Arial';
    c.FontSize  = 12;
    ax          = gca;
    ax.FontSize = 16;
    ax.LineWidth = 2;
    ax.Position = position(ii,:);
    X_Tick     = 5:12:81;
    ax.XTick   = X_Tick;
    % Map longitude labels
    X_Tick_Label_Pre = [];
    for i = I6:I7
        X_Tick_Label_Pre{i-I6+1} = Lon_Label{i};
    end
    X_Tick_Label = [];
    for i = 1:length(X_Tick)
        ax.XTickLabel(i) = X_Tick_Label_Pre{X_Tick(i)};
    end
    ax.YTick       = [1 11 21 31 41 51 61];
    ax.YTickLabel  = {'-30','-20','-10','0','10','20','30'};
    ax.FontName    = 'Arial';

    % Optional: grid settings (commented out)
    %grid on
    %set(gca,'GridLineStyle','--','GridColor','k','GridAlpha',.4);
    %ax.GridLineWid = .4;
end
