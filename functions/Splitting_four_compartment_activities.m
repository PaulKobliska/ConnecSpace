function      [spkx_in1,spky_in1,spkx_in2,spky_in2,spkx_in3,spky_in3,spkx_in4,spky_in4,spkx,spky,ts_in1,ts_in2,ts_in3,ts_in4,posgx_in1,posgy_in1,posgx_in2,posgy_in2,posgx_in3,posgy_in3,posgx_in4,posgy_in4,rate_in1,rate_in2,rate_in3,rate_in4,ratemap_in1,ratemap_in2,ratemap_in3,ratemap_in4,rotated_ratemap_in1,rotated_ratemap_in2,rotated_ratemap_in3,rotated_ratemap_in4,RrankO1,RrankO2,RrankO3,RrankO4,RrankO5,RrankO6,rotated_RrankO1,rotated_RrankO2,rotated_RrankO3,rotated_RrankO4,rotated_RrankO5,rotated_RrankO6] = Splitting_four_compartment_activities(ts,posts,posrx,posry,posgx,posgy,posrx_in1,posry_in1,posrx_in2,posry_in2,posrx_in3,posry_in3,posrx_in4,posry_in4,x,x2,x3,x4,y,y2,y3,y4,pos_in_rectangle1,pos_in_rectangle2,pos_in_rectangle3,pos_in_rectangle4,posts_in1,posts_in2,posts_in3,posts_in4,binWidth,shape,V,V2,V3,V4,xx,xx2,xx3,xx4);


N = length(ts);
spkx = zeros(N,1);
spky = zeros(N,1);
for i = 1:N
    [val,ind] = min(abs(posts-ts(i)));
    spkx(i) = posrx(ind); 
    spky(i) = posry(ind);
end

%% Isolate spike position in the 1st compartment (up left)
spike_in_rectangle1_x=zeros(length(spkx),1);spike_in_rectangle1_x=logical(spike_in_rectangle1_x);
spike_in_rectangle1_y=zeros(length(spkx),1);spike_in_rectangle1_y=logical(spike_in_rectangle1_y);

for p = 1:length(V)
    a=spkx<xx(p) & spky>V(p);
    spike_in_rectangle1_x(a) = 1;
end

for q = 1:length(V2)
    a2=spkx<xx2(q) & spky>V2(q);
    spike_in_rectangle1_y(a2) = 1;
end
spike_in_rectangle1(find(spike_in_rectangle1_x==1 & spike_in_rectangle1_y==1))=1;

spike_in_rectangle1=logical(spike_in_rectangle1);
spkx_in1 = spkx(spike_in_rectangle1);
spky_in1 = spky(spike_in_rectangle1);
ts_in1 = ts(spike_in_rectangle1);
rate_in1 =(length(spkx_in1)/(length(posrx_in1)*40))*1000;
posgx_in1 = posgx(pos_in_rectangle1);
posgy_in1 = posgy(pos_in_rectangle1);

% scatter(spkx_in1,spky_in1);
% scatter(spkx,spky,'x');

%% plot rate map in rectangle symetrie axiale
[posrx_in1,posry_in1] = smooth_path(posrx_in1,posry_in1);
[posrx_in1,posry_in1] = center_path(posrx_in1,posry_in1,shape);

[mapAxis_in1] = mapaxis(posrx_in1,posry_in1,binWidth);
visited_in1 = visitedBins(posrx_in1,posry_in1,mapAxis_in1); % calculate the visited zones of the arena 
[spkx_in1,spky_in1] = get_pos_spikes(ts_in1,posrx_in1,posry_in1,posts_in1); % get the position of the spikes

[ratemap_in1] = ratemap_gaussian(2*binWidth,spkx_in1,spky_in1,posrx_in1,posry_in1,posts_in1,binWidth,mapAxis_in1); % build matrix ratemap 
% remove invisited bins from the rate map
ratemap_in1(visited_in1==0) = NaN;


%% plot rate map in rectangle symétrie radiale
minus_posrx_in1 = -(posrx_in1);
minus_posry_in1 = -(posry_in1);

[minus_posrx_in1,minus_posry_in1] = smooth_path(minus_posrx_in1,minus_posry_in1);
[minus_posrx_in1,minus_posry_in1] = center_path(minus_posrx_in1,minus_posry_in1,shape);

[rotated_mapAxis_in1] = mapAxis_in1;%mapaxis(minus_posrx_in,minus_posry_in,binWidth);
rotated_visited_in1 = visitedBins(minus_posrx_in1,minus_posry_in1,rotated_mapAxis_in1); % calculate the visited zones of the arena 
[spkx_in1,spky_in1] = get_pos_spikes(ts_in1,minus_posrx_in1,minus_posry_in1,posts_in1); % get the position of the spikes

[rotated_ratemap_in1] = ratemap_gaussian(2*binWidth,spkx_in1,spky_in1,minus_posrx_in1,minus_posry_in1,posts_in1,binWidth,rotated_mapAxis_in1); % build matrix ratemap 
% remove invisited bins from the rate map
rotated_ratemap_in1(rotated_visited_in1==0) = NaN;

%% Isolate spike position in the 2nd compartment (down left)
spike_in_rectangle2_x=zeros(length(spkx),1);spike_in_rectangle2_x=logical(spike_in_rectangle2_x);
spike_in_rectangle2_y=zeros(length(spkx),1);spike_in_rectangle2_y=logical(spike_in_rectangle2_y);

for p = 1:length(V3)
    b=spkx<xx3(p) & spky<=V3(p);
    spike_in_rectangle2_x(b) = 1;
end

for q = 1:length(V2)
    b2=spkx<=xx2(q) & spky<V2(q);
    spike_in_rectangle2_y(b2) = 1;
end
spike_in_rectangle2(find(spike_in_rectangle2_x==1 & spike_in_rectangle2_y==1))=1;

spike_in_rectangle2=logical(spike_in_rectangle2);
spkx_in2 = spkx(spike_in_rectangle2);
spky_in2 = spky(spike_in_rectangle2);
ts_in2 = ts(spike_in_rectangle2);
rate_in2 =(length(spkx_in2)/(length(posrx_in2)*40))*1000;
posgx_in2 = posgx(pos_in_rectangle2);
posgy_in2 = posgy(pos_in_rectangle2);

% scatter(spkx_in2,spky_in2);
% scatter(spkx,spky,'x');

%% plot rate map in rectangle symetrie axiale
[posrx_in2,posry_in2] = smooth_path(posrx_in2,posry_in2);
[posrx_in2,posry_in2] = center_path(posrx_in2,posry_in2,shape);

[mapAxis_in2] = mapaxis(posrx_in2,posry_in2,binWidth);
visited_in2 = visitedBins(posrx_in2,posry_in2,mapAxis_in2); % calculate the visited zones of the arena 
[spkx_in2,spky_in2] = get_pos_spikes(ts_in2,posrx_in2,posry_in2,posts_in2); % get the position of the spikes

[ratemap_in2] = ratemap_gaussian(2*binWidth,spkx_in2,spky_in2,posrx_in2,posry_in2,posts_in2,binWidth,mapAxis_in2); % build matrix ratemap 
% remove invisited bins from the rate map
ratemap_in2(visited_in2==0) = NaN;


%% plot rate map in rectangle symétrie radiale
minus_posrx_in2 = -(posrx_in2);
minus_posry_in2 = -(posry_in2);

[minus_posrx_in2,minus_posry_in2] = smooth_path(minus_posrx_in2,minus_posry_in2);
[minus_posrx_in2,minus_posry_in2] = center_path(minus_posrx_in2,minus_posry_in2,shape);

[rotated_mapAxis_in2] = mapAxis_in2;%mapaxis(minus_posrx_in,minus_posry_in,binWidth);
rotated_visited_in2 = visitedBins(minus_posrx_in2,minus_posry_in2,rotated_mapAxis_in2); % calculate the visited zones of the arena 
[spkx_in2,spky_in2] = get_pos_spikes(ts_in2,minus_posrx_in2,minus_posry_in2,posts_in2); % get the position of the spikes

[rotated_ratemap_in2] = ratemap_gaussian(2*binWidth,spkx_in2,spky_in2,minus_posrx_in2,minus_posry_in2,posts_in2,binWidth,rotated_mapAxis_in2); % build matrix ratemap 
% remove invisited bins from the rate map
rotated_ratemap_in2(rotated_visited_in2==0) = NaN;

%% Isolate spike position in the 3rd compartment (down right)
spike_in_rectangle3_x=zeros(length(spkx),1);spike_in_rectangle3_x=logical(spike_in_rectangle3_x);
spike_in_rectangle3_y=zeros(length(spkx),1);spike_in_rectangle3_y=logical(spike_in_rectangle3_y);

for p = 1:length(V4)
    c=spkx>=xx4(p) & spky<V4(p);
    spike_in_rectangle3_x(c) = 1;
end

for q = 1:length(V3)
    c2=spkx>xx3(q) & spky<=V3(q);
    spike_in_rectangle3_y(c2) = 1;
end
spike_in_rectangle3(find(spike_in_rectangle3_x==1 & spike_in_rectangle3_y==1))=1;

spike_in_rectangle3=logical(spike_in_rectangle3);
spkx_in3 = spkx(spike_in_rectangle3);
spky_in3 = spky(spike_in_rectangle3);
ts_in3 = ts(spike_in_rectangle3);
rate_in3 =(length(spkx_in3)/(length(posrx_in3)*40))*1000;
posgx_in3 = posgx(pos_in_rectangle3);
posgy_in3 = posgy(pos_in_rectangle3);

% scatter(spkx_in3,spky_in3);
% scatter(spkx,spky,'x');

%% plot rate map in rectangle symetrie axiale
[posrx_in3,posry_in3] = smooth_path(posrx_in3,posry_in3);
[posrx_in3,posry_in3] = center_path(posrx_in3,posry_in3,shape);

[mapAxis_in3] = mapaxis(posrx_in3,posry_in3,binWidth);
visited_in3 = visitedBins(posrx_in3,posry_in3,mapAxis_in3); % calculate the visited zones of the arena 
[spkx_in3,spky_in3] = get_pos_spikes(ts_in3,posrx_in3,posry_in3,posts_in3); % get the position of the spikes

[ratemap_in3] = ratemap_gaussian(2*binWidth,spkx_in3,spky_in3,posrx_in3,posry_in3,posts_in3,binWidth,mapAxis_in3); % build matrix ratemap 
% remove invisited bins from the rate map
ratemap_in3(visited_in3==0) = NaN;


%% plot rate map in rectangle symétrie radiale
minus_posrx_in3 = -(posrx_in3);
minus_posry_in3 = -(posry_in3);

[minus_posrx_in3,minus_posry_in3] = smooth_path(minus_posrx_in3,minus_posry_in3);
[minus_posrx_in3,minus_posry_in3] = center_path(minus_posrx_in3,minus_posry_in3,shape);

[rotated_mapAxis_in3] = mapAxis_in3;%mapaxis(minus_posrx_in,minus_posry_in,binWidth);
rotated_visited_in3 = visitedBins(minus_posrx_in3,minus_posry_in3,rotated_mapAxis_in3); % calculate the visited zones of the arena 
[spkx_in3,spky_in3] = get_pos_spikes(ts_in3,minus_posrx_in3,minus_posry_in3,posts_in3); % get the position of the spikes

[rotated_ratemap_in3] = ratemap_gaussian(2*binWidth,spkx_in3,spky_in3,minus_posrx_in3,minus_posry_in3,posts_in3,binWidth,rotated_mapAxis_in3); % build matrix ratemap 
% remove invisited bins from the rate map
rotated_ratemap_in3(rotated_visited_in3==0) = NaN;

%% Isolate spike position in the 4th compartment (down right)
spike_in_rectangle4_x=zeros(length(spkx),1);spike_in_rectangle4_x=logical(spike_in_rectangle4_x);
spike_in_rectangle4_y=zeros(length(spkx),1);spike_in_rectangle4_y=logical(spike_in_rectangle4_y);

for p = 1:length(V)
    d=spkx>xx(p) & spky>V(p);
    spike_in_rectangle4_x(d) = 1;
end

for q = 1:length(V4)
    d2=spkx>=xx4(q) & spky>V4(q);
    spike_in_rectangle4_y(d2) = 1;
end
spike_in_rectangle4(find(spike_in_rectangle4_x==1 & spike_in_rectangle4_y==1))=1;

spike_in_rectangle4=logical(spike_in_rectangle4);
spkx_in4 = spkx(spike_in_rectangle4);
spky_in4 = spky(spike_in_rectangle4);
ts_in4 = ts(spike_in_rectangle4);
rate_in4 =(length(spkx_in4)/(length(posrx_in4)*40))*1000;
posgx_in4 = posgx(pos_in_rectangle4);
posgy_in4 = posgy(pos_in_rectangle4);

% scatter(spkx_in4,spky_in4);
% scatter(spkx,spky,'x');

%% plot rate map in rectangle symetrie axiale
[posrx_in4,posry_in4] = smooth_path(posrx_in4,posry_in4);
[posrx_in4,posry_in4] = center_path(posrx_in4,posry_in4,shape);

[mapAxis_in4] = mapaxis(posrx_in4,posry_in4,binWidth);
visited_in4 = visitedBins(posrx_in4,posry_in4,mapAxis_in4); % calculate the visited zones of the arena 
[spkx_in4,spky_in4] = get_pos_spikes(ts_in4,posrx_in4,posry_in4,posts_in4); % get the position of the spikes

[ratemap_in4] = ratemap_gaussian(2*binWidth,spkx_in4,spky_in4,posrx_in4,posry_in4,posts_in4,binWidth,mapAxis_in4); % build matrix ratemap 
% remove invisited bins from the rate map
ratemap_in4(visited_in4==0) = NaN;


%% plot rate map in rectangle symétrie radiale
minus_posrx_in4 = -(posrx_in4);
minus_posry_in4 = -(posry_in4);

[minus_posrx_in4,minus_posry_in4] = smooth_path(minus_posrx_in4,minus_posry_in4);
[minus_posrx_in4,minus_posry_in4] = center_path(minus_posrx_in4,minus_posry_in4,shape);

[rotated_mapAxis_in4] = mapAxis_in4;%mapaxis(minus_posrx_in,minus_posry_in,binWidth);
rotated_visited_in4 = visitedBins(minus_posrx_in4,minus_posry_in4,rotated_mapAxis_in4); % calculate the visited zones of the arena 
[spkx_in4,spky_in4] = get_pos_spikes(ts_in4,minus_posrx_in4,minus_posry_in4,posts_in4); % get the position of the spikes

[rotated_ratemap_in4] = ratemap_gaussian(2*binWidth,spkx_in4,spky_in4,minus_posrx_in4,minus_posry_in4,posts_in4,binWidth,rotated_mapAxis_in4); % build matrix ratemap 
% remove invisited bins from the rate map
rotated_ratemap_in4(rotated_visited_in4==0) = NaN;
%% Compartment correlation
minLength = min([length(ratemap_in1), length(ratemap_in2), length(ratemap_in3), length(ratemap_in4)]);
ratemap_in1=ratemap_in1(1:minLength,1:minLength); ratemap_in2=ratemap_in2(1:minLength,1:minLength); ratemap_in3=ratemap_in3(1:minLength,1:minLength); ratemap_in4=ratemap_in4(1:minLength,1:minLength); 
 
[corrcoeff1,nonZero_corrcoeff1,RrankO1] = get_2D_correlation_value(ratemap_in1,ratemap_in2);
[corrcoeff2,nonZero_corrcoeff2,RrankO2] = get_2D_correlation_value(ratemap_in1,ratemap_in3);
[corrcoeff3,nonZero_corrcoeff3,RrankO3] = get_2D_correlation_value(ratemap_in1,ratemap_in4);
[corrcoeff4,nonZero_corrcoeff4,RrankO4] = get_2D_correlation_value(ratemap_in2,ratemap_in3);
[corrcoeff5,nonZero_corrcoeff5,RrankO5] = get_2D_correlation_value(ratemap_in2,ratemap_in4);
[corrcoeff6,nonZero_corrcoeff6,RrankO6] = get_2D_correlation_value(ratemap_in3,ratemap_in4);

minLength = min([length(rotated_ratemap_in1), length(rotated_ratemap_in2), length(rotated_ratemap_in3), length(rotated_ratemap_in4)]);
rotated_ratemap_in1=rotated_ratemap_in1(1:minLength,1:minLength); rotated_ratemap_in2=rotated_ratemap_in2(1:minLength,1:minLength); rotated_ratemap_in3=rotated_ratemap_in3(1:minLength,1:minLength); rotated_ratemap_in4=rotated_ratemap_in4(1:minLength,1:minLength); 
[rotated_corrcoeff1,rotated_nonZero_corrcoeff1,rotated_RrankO1] = get_2D_correlation_value(rotated_ratemap_in1,rotated_ratemap_in2);
[rotated_corrcoeff2,rotated_nonZero_corrcoeff2,rotated_RrankO2] = get_2D_correlation_value(rotated_ratemap_in1,rotated_ratemap_in3);
[rotated_corrcoeff3,rotated_nonZero_corrcoeff3,rotated_RrankO3] = get_2D_correlation_value(rotated_ratemap_in1,rotated_ratemap_in4);
[rotated_corrcoeff4,rotated_nonZero_corrcoeff4,rotated_RrankO4] = get_2D_correlation_value(rotated_ratemap_in2,rotated_ratemap_in3);
[rotated_corrcoeff5,rotated_nonZero_corrcoeff5,rotated_RrankO5] = get_2D_correlation_value(rotated_ratemap_in2,rotated_ratemap_in4);
[rotated_corrcoeff6,rotated_nonZero_corrcoeff6,rotated_RrankO6] = get_2D_correlation_value(rotated_ratemap_in3,rotated_ratemap_in4);
