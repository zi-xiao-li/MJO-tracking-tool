% Track MJO phase speed within the slices

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

% Find the index corresponding to the year 
Number_All = [1979:2013];
Number     = find(Number_All==The_Chosen_One);

ex_OLR     = Segment_OLR(:,:,Number)';
ex_x0      = 90;   % Corresponding longitude  
ex_time    = Segment_Time(:,[1:3],Number);
ex_lon     = Segment_Lon;
EX_t0      = Segment_t0{Number}; % That is, t0 for 2006, in year-month-day format

% If no MJO events occurred in that year >>>>>>>>>>>>>>>>>>>>>>>>
if isempty(EX_t0)==1
    continue
end

















for H=1:size(EX_t0,1)



    ex_t0(1,:) = EX_t0(H,:);

% In the 181-day Hovmöller plot, y-axis: t ∈ [t0-12 : t0+12], coordinate of t0
ex_ref_y0  = find(ex_time(:,1)==ex_t0(1,1) & ex_time(:,2)==ex_t0(1,2) & ex_time(:,3)==ex_t0(1,3));
% Skip out-of-bounds events, e.g., May 29
% *This is reasonable because we selected dates from October, November, December, January, February, March, April, May; too distant dates should not matter
% *At worst, we can remove them later
%! Previous out-of-bounds issues occurred because Segment_t0 was selected in May and October; after deletion, this problem does not appear
%! The slice selection in May and October is itself intended to preempt out-of-bounds problems
% if ex_ref_y0 - 12 < 1 | ex_ref_y0 + 12 > size(ex_OLR, 1)
%     disp(['Skipping out-of-bounds event: ', num2str(ex_t0(1,:))])
%     continue
% end
ex_y0      = [ex_ref_y0-12:ex_ref_y0+12]-0.5; % 12 days before and after the reference y0



% If t ∈ [t0-12, t0+12] contains values less than 0, e.g., reference point on Nov 1 >>>>>>>>>>>>>>>>>>>>>>>>
% if isempty(find(ex_y0 < 0)) == 0
%     ex_y0 = ex_y0(find(ex_y0 > 0)); % Remove 0 values because later we use y0-0.5, x0-0.5; this gives the y-axis coordinates
% end



% If t ∈ [t0-12, t0+12] contains values greater than 181, e.g., reference point on Apr 30 >>>>>>>>>>>>>>>>>>>>>>>>
% if isempty(find(ex_y0 > 181)) == 0
%     ex_y0 = ex_y0(find(ex_y0 < 181)); % Remove values exceeding 181 because later we use y0-0.5, x0-0.5; this gives the y-axis coordinates
% end




% Find the line using t ∈ [y0-12, y0+12]
syms x y

[~, x0]    = min(abs(ex_lon(:) - ex_x0)); % The reference point longitude is fixed, i.e., 90°E
x0          = x0 - 0.5;
% Calculate K values, every 0.1 m/s
K_real = 1:0.1:25;
K      = 111001 ./ (34560 .* K_real + 1);


% From here, start displaying on the coordinate axes >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

% Determine the intersection points of each line on the right half with the grid
coords_left  = [];
coords_right = [];

coords_right = zeros(length(ex_y0), length(K), 2); % For each slope passing through the reference point within the date range, find the intersection with the grid edge
for ii = 1:length(ex_y0)
    y0 = ex_y0(ii);
    for jj = 1:length(K)
        k  = K(jj);
        y = @(x) k*(x - x0) + y0;
        if (y(81) <= 243) % If the slope does not exceed the grid, coordinates are (81, y(81))
            coords_right(ii,jj,1) = 81;
            coords_right(ii,jj,2) = y(81);
        else
            y = @(x) (1./k)*(x - y0) + x0; % If it exceeds, set coordinates to (x(181), 181)
            coords_right(ii,jj,1) = y(243);
            coords_right(ii,jj,2) = 243;
        end
    end
end

% Determine the intersection points of each line on the left half with the grid
coords_left = zeros(length(ex_y0), length(K), 2); 
for ii = 1:length(ex_y0)
    y0 = ex_y0(ii);
    for jj = 1:length(K)
        k  = K(jj);
        y = @(x) k*(x - x0) + y0;
        if (y(0) >= 0) % If the slope does not exceed the grid, coordinates are (0, y(0))
            coords_left(ii,jj,1) = 0;
            coords_left(ii,jj,2) = y(0);
        else 
            y = @(x) (1./k)*(x - y0) + x0; % If it exceeds, set coordinates to (x(0), 0)
            coords_left(ii,jj,1) = y(0);
            coords_left(ii,jj,2) = 0;
        end
    end
end

% Calculate the grid cells that the line passes through
% For a given day and slope, with known left and right vertices, determine the region the line crosses
segs_info = [];
segs_info = cell(size(coords_left,1), size(coords_left,2));

for ii = 1:size(coords_right,1)
    y0 = ex_y0(ii);
    for jj = 1:size(coords_right,2)
        k     = [];
        P1    = [];
        P2    = [];
        xdx   = [];
        xdy   = [];
        ydx   = [];
        ydy   = [];
        k        = K(jj);
        
        % Identify the confirmed start and end points
        P1(:,:)  = coords_left(ii,jj,:);
        P2(:,:)  = coords_right(ii,jj,:);
        P1_x     = P1(1);
        P1_y     = P1(2);
        P2_x     = P2(1);
        P2_y     = P2(2);
        xmin     = min(P1_x, P2_x);
        ymin     = min(P1_y, P2_y);
        xmax     = max(P1_x, P2_x);
        ymax     = max(P1_y, P2_y);
        
        % Round x to find corresponding y coordinates
        y = @(x) k*(x - x0) + y0;
        for i = ceil(xmin):floor(xmax)
            xdx(i - ceil(xmin) + 1) = i;
            xdy(i - ceil(xmin) + 1) = y(i);
            if (y(i) < 0 && y(i) > -1)
                xdy(i - ceil(xmin) + 1) = 0;
                disp(['Value below 0: ', num2str(y(i)), ' Reference point: ', num2str(x0), ',', num2str(y0), ...
                        ' Day index: ', num2str(ii), ', Slope index: ', num2str(jj), ', XY coords: ', num2str(i), ',', num2str(y(i))]);
            elseif (y(i) < -1)
                disp(['Value below -1: ', num2str(y(i)), ' Reference point: ', num2str(x0), ',', num2str(y0), ...
                        ' Day index: ', num2str(ii), ', Slope index: ', num2str(jj), ', XY coords: ', num2str(i), ',', num2str(y(i))]);
            end
        end

        % Round y to find corresponding x coordinates
        y = @(x) (1./k)*(x - y0) + x0;
        for j = ceil(ymin):floor(ymax)
            ydy(j - ceil(ymin) + 1) = j;
            ydx(j - ceil(ymin) + 1) = y(j);
            if (y(j) < 0 && y(j) > -1)
                ydx(j - ceil(ymin) + 1) = 0;
                disp(['Value below 0: ', num2str(y(j)), ' Reference point: ', num2str(x0), ',', num2str(y0), ...
                        ' Day index: ', num2str(ii), ', Slope index: ', num2str(jj), ', XY coords: ', num2str(y(j)), ',', num2str(j)]);
            elseif (y(j) < -1)
                disp(['Value below -1: ', num2str(y(j)), ' Reference point: ', num2str(x0), ',', num2str(y0), ...
                        ' Day index: ', num2str(ii), ', Slope index: ', num2str(jj), ', XY coords: ', num2str(y(j)), ',', num2str(j)]);
            end
        end

        % Merge all coordinates and remove duplicates
        SP = unique([P1_x, P2_x, xdx, ydx; P1_y, P2_y, xdy, ydy]','rows'); % Monotonically increasing function

        segs = [];
        for t = 1 : size(SP,1)-1
            segs(t).index_x = max(ceil(SP(t+1,1)), ceil(SP(t,1))); % Column index
            segs(t).index_y = max(ceil(SP(t+1,2)), ceil(SP(t,2))); % Row index
            segs(t).OLRA    = ex_OLR(segs(t).index_y, segs(t).index_x);
        end
        
        % Store the line segments passing through grid points
        segs_info{ii,jj} = segs;
    end
end

% Identify the longest segment along each trail line where OLRA is below one standard deviation under the corresponding mean.
% Segments with longitude gaps < 10° are considered continuous events.
% Only integer longitudes are considered.

for ii = 1:size(segs_info,1)
    for jj = 1:size(segs_info,2)
        var1   = [];
        Series = [];
        var1   = segs_info{ii,jj};

        % Identify all grid points with OLRA below mean minus 1 SD
        for t = 1:length(var1)
            if var1(t).OLRA <= Segment_mean(find(Segment_std(:,2) == var1(t).index_x)) ...
                                - Segment_std(find(Segment_std(:,2) == var1(t).index_x))
                var1(t).track = 1;
            else 
                var1(t).track = 0;
            end
            Series(t) = var1(t).track;
        end

        % Merge zero segments separated by less than 10° longitude
        temp = []; aa_start = []; aa_end = [];
        temp = diff(Series);
        aa_start = find(temp == -1) + 1; % Start of zero segments
        aa_end   = find(temp == 1);      % End of zero segments

        if Series(1) == 0
            aa_start = [1 aa_start];
        end
        if Series(end) == 0
            aa_end = [aa_end length(Series)];
        end

        for i = 1:length(aa_end)
            if var1(aa_end(i)).index_x - var1(aa_start(i)).index_x <= 4
                Series(aa_start(i):aa_end(i)) = 1;
            end
        end

        % Update tracking for all points
        for t = 1:length(segs_info{ii,jj})
            segs_info{ii,jj}(t).track = Series(t);
        end

        % Find the longest sequence of 1s (longest segment), ignoring small fragments
        temp = []; aa_start = []; aa_end = [];
        temp = diff(Series);
        aa_start = find(temp == 1) + 1; % Start of 1 segments
        aa_end   = find(temp == -1);    % End of 1 segments

        if Series(1) == 1
            aa_start = [1 aa_start];
        end
        if Series(end) == 1
            aa_end = [aa_end length(Series)];
        end

        % Calculate A(t,c) and L(t,c)
        a = []; b = []; bb_start = []; bb_end = []; cc = 0;
        if isempty(aa_end)
            segs_info{ii,jj}(1).StartLong = 0;
            segs_info{ii,jj}(1).EndLong   = 0;
            segs_info{ii,jj}(1).A         = 0;
            segs_info{ii,jj}(1).L         = 0;
            segs_info{ii,jj}(1).Start     = 0;
            segs_info{ii,jj}(1).End       = 0;
        else
            [a, b] = max(aa_end - aa_start);
            bb_start = aa_start(b);
            bb_end   = aa_end(b); 
            for t = bb_start:bb_end
                 cc = cc + segs_info{ii,jj}(t).OLRA;
            end
            segs_info{ii,jj}(1).StartLong = segs_info{ii,jj}(bb_start).index_x;
            segs_info{ii,jj}(1).EndLong   = segs_info{ii,jj}(bb_end).index_x;
            segs_info{ii,jj}(1).StartDate = segs_info{ii,jj}(bb_start).index_y;
            segs_info{ii,jj}(1).EndDate   = segs_info{ii,jj}(bb_end).index_y;
            segs_info{ii,jj}(1).A         = cc;
            segs_info{ii,jj}(1).L         = 2.5 * (segs_info{ii,jj}(bb_end).index_x - segs_info{ii,jj}(bb_start).index_x);
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

% Find A_m and L_m
A_m = max(max(A));
L_m = max(max(L));

% Calculate B(t,c)
B = A ./ A_m + L ./ L_m;
B_m = max(max(B));
[pos_x, pos_y] = find(B == B_m);
% Corresponding date and slope: pos_x for date, pos_y for slope

% If K slope spacing is too small, multiple pos_y values may appear, same for dates >>>>>>>>>>>>>>>>>>>>>>>>

if length(pos_x) ~= 1
    disp(['pos x (date) ', num2str(pos_x')])
    pos_x = max(pos_x); % May result in a longer lifecycle
    % pos_x = nan;
end

if length(pos_y) ~= 1
    disp(['pos y (slope) ', num2str(pos_y')])
    pos_y = max(pos_y);  % Slower speed
    % pos_y = nan;
end

Final_K    = [];
Final_info = [];
Final_K    = K_real(pos_y);
Final_info = segs_info{pos_x,pos_y};


disp(['Start: ', num2str(ex_time(Final_info(1).StartDate,:))])
disp(['End: ', num2str(ex_time(Final_info(1).EndDate,:))])
disp(['Day 0: ', num2str(ex_t0)])
disp(['Start Longitude: ', num2str(Final_info(1).StartLong * 2.5 + 20)])
disp(['End Longitude: ', num2str(Final_info(1).EndLong * 2.5 + 20)])
disp(['MJO propagation speed: ', num2str(Final_K), ' m/s'])

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




file_name1 = "D:\project\output\PropagationInfo.xlsx";  %% Note: remove duplicate values and incorrect dates
file_name2 = "D:\project\output\PropagationInfo_prepare.xlsx";

xlswrite(file_name1, data_print)

% Read data
data = [];
data = xlsread(file_name1);

% Remove certain values from the data
aa1  = [];
aa1  = find(data(:,11) <= 80 & data(:,12) >= 120);
data = data(aa1, :);

xlswrite(file_name2, data)

