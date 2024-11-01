function  [xy_left,xy_right,xy_door,left,right,door] = Splitting_compartment_v1_removing_door(PosMtx,inputfile)
%modified by PYves 02/09/2022 : draw polygon of compartments and door.
%coordonnees pour cellules directionnelles
shape = 'circle';
posgx = PosMtx(:,2);
posgy = PosMtx(:,3);
posrx = PosMtx(:,4);
posry = PosMtx(:,5);

posts = PosMtx(:,1);
% correct the x y coordinates 
%[posgx,posgy,posts] = remJumpTrack(posgx,posgy,posts,threshold);
% [posgx,posgy,posts] = interpolatePos(posgx,posgy,posts);
% then smooth and center
[posgx,posgy] = smooth_path(posgx,posgy);
[posgx,posgy] = center_path(posgx,posgy,shape);

posts = PosMtx(:,1);
% posts = OpenfIeld.pos.ts;
% correct the x y coordinates
%[posrx,posry,posts] = remJumpTrack(posrx,posry,posts,threshold);
% [posrx,posry,posts] = interpolatePos(posrx,posry,posts);
% then smooth and center
[posrx,posry] = smooth_path(posrx,posry);
[posrx,posry] = center_path(posrx,posry,shape);

figure
plot(posrx,posry,'color',[.6 .6 .6]); set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
hold on
plot(posgx,posgy,'color',[.2 .2 .2]);

roi_left = 	drawpolygon;
roi_right = drawpolygon;
roi_door = 	drawpolygon;

xy_left = roi_left.Position;
xy_right = roi_right.Position;
xy_door = roi_door.Position;

[left,~] = inpolygon(posrx,posry,roi_left.Position(:,1),roi_left.Position(:,2));
[right,~] = inpolygon(posrx,posry,roi_right.Position(:,1),roi_right.Position(:,2));
[door,~] = inpolygon(posrx,posry,roi_door.Position(:,1),roi_door.Position(:,2));

plot(posrx,posry); set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
hold on
% plot position of position in each quadrant with different color
scatter(posrx(left),posry(left),'r.'); % points inside
scatter(posrx(right),posry(right),'g.'); % points inside
scatter(posrx(door),posry(door),'b.'); 


% % % % % % % [x,y] = ginput(2); sep = line([x(1) x(2)], [y(1) y(2)],'Color','r');
% % % % % % % % hold on
% % % % % % % % scatter(x,y,'filled');
% % % % % % % 
% % % % % % % xx=min(x):0.1:max(x); 
% % % % % % % coefficients=polyfit([x(1), x(2)],[y(1),y(2)],1);
% % % % % % % V=coefficients(1)*xx+coefficients(2);
% % % % % % % % plot(xx,V)
% % % % % % % 
% % % % % % % minx = min(posrx)-30; miny = min(posry)-30; maxx = max(posrx)+30; maxy = max(posry)+30;
% % % % % % % lengthx = max(x) - minx;
% % % % % % % lengthy = max(y) - miny;
% % % % % % % % hold on
% % % % % % % % rectangle('Position',[minx miny lengthx lengthy] )
% % % % % % % V=V'; xx=xx';
% % % % % % % pos_in_rectangle=zeros(length(posrx),1);
% % % % % % % pos_in_rectangle=logical(pos_in_rectangle);
% % % % % % % for j = 1:length(V)
% % % % % % % a = ((posrx<min(xx(j))) & (posry<V(j)));
% % % % % % % pos_in_rectangle(a)=1;
% % % % % % % end
% % % % % % % pos_out_rectangle = pos_in_rectangle == false;
% % % % % % % 
% % % % % % % % eliminate coordinate zero from position
% % % % % % % pos_in_rectangle(PosMtx(:,2)==0,:) = []; 
% % % % % % % pos_in_rectangle(PosMtx(:,3)==0,:) = [];
% % % % % % % pos_out_rectangle(PosMtx(:,2)==0,:) = []; 
% % % % % % % pos_out_rectangle(PosMtx(:,3)==0,:) = [];
% % % % % % % 
% % % % % % % posrx_in = posrx(pos_in_rectangle);
% % % % % % % posry_in = posry(pos_in_rectangle);
% % % % % % % posts_in = posts(pos_in_rectangle);
% % % % % % % posrx_out = posrx(pos_out_rectangle);
% % % % % % % posry_out = posry(pos_out_rectangle);
% % % % % % % posts_out = posts(pos_out_rectangle);
% % % % % % % 
% % % % % % % hold on
% % % % % % % plot(posrx_in,posry_in);
% % % % % % % plot(posrx_out,posry_out);
% % % % % % % 
% % % % % % % 
% % % % % % % hold on
% % % % % % % [posdoory,posdoorx] = ginput(1);
end

