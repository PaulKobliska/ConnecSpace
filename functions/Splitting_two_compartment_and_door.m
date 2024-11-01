function  [spkx,spky,spike_left,spike_right,spike_door] = Splitting_two_compartment_and_door(ts,posx,posy,posts,xy_left,xy_right,xy_door)

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
[spike_left,~] = inpolygon(spkx,spky,xy_left(:,1),xy_left(:,2));
[spike_right,~] = inpolygon(spkx,spky,xy_right(:,1),xy_right(:,2));
[spike_door,~] = inpolygon(spkx,spky,xy_door(:,1),xy_door(:,2));




 
