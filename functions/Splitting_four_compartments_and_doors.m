function  [spkx,spky,spike_top_left,spike_bottom_left,spike_top_right,spike_bottom_right,door_top,door_left,door_bottom,door_right] = Splitting_four_compartments_and_doors(ts,posx,posy,posts,xy_top_left,xy_bottom_left,xy_top_right,xy_bottom_right,xy_top_door,xy_left_door,xy_bottom_door,xy_right_door)
%% spike extrac for each quadrant
N = length(ts);
spkx = zeros(N,1);
spky = zeros(N,1);
for i = 1:N
    [~,ind] = min(abs(posts-ts(i)));
    spkx(i) = posx(ind);
    spky(i) = posy(ind);
end
% extract spikes coordinates of each quadrant
[spike_top_left,~] = inpolygon(spkx,spky,xy_top_left(:,1),xy_top_left(:,2));
[spike_bottom_left,~] = inpolygon(spkx,spky,xy_bottom_left(:,1),xy_bottom_left(:,2));
[spike_top_right,~] = inpolygon(spkx,spky,xy_top_right(:,1),xy_top_right(:,2));
[spike_bottom_right,~] = inpolygon(spkx,spky,xy_bottom_right(:,1),xy_bottom_right(:,2));

[door_top,~] = inpolygon(spkx,spky,xy_top_door(:,1),xy_top_door(:,2));
[door_left,~] = inpolygon(spkx,spky,xy_left_door(:,1),xy_left_door(:,2));
[door_bottom,~] = inpolygon(spkx,spky,xy_bottom_door(:,1),xy_bottom_door(:,2));
[door_right,~] = inpolygon(spkx,spky,xy_right_door(:,1),xy_right_door(:,2));




