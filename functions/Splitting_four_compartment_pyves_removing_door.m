function  [xy_top_left,xy_bottom_left,xy_top_right,xy_bottom_right,xy_top_door,xy_left_door,xy_bottom_door,xy_right_door,topleft,bottomleft,bottomright,topright,topdoor,leftdoor,bottomdoor,rightdoor] = Splitting_four_compartment_pyves_removing_door(PosMtx);
%modified by PYves 02/09/2022 : draw polygon of compartments and door.

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

roi_top_left = 	drawpolygon;
roi_bottom_left = drawpolygon;
roi_bottom_right = 	drawpolygon;
roi_top_right = 	drawpolygon;

roi_top_door = 	drawpolygon;
roi_left_door = 	drawpolygon;
roi_bottom_door = 	drawpolygon;
roi_right_door = 	drawpolygon;

xy_top_left = roi_top_left.Position;
xy_bottom_left = roi_bottom_left.Position;
xy_top_right = roi_top_right.Position;
xy_bottom_right = roi_bottom_right.Position;

xy_top_door = roi_top_door.Position;
xy_left_door = roi_left_door.Position;
xy_bottom_door = roi_bottom_door.Position;
xy_right_door = roi_right_door.Position;



% rect = getrect;

[topleft,~] = inpolygon(posrx,posry,roi_top_left.Position(:,1),roi_top_left.Position(:,2));
[bottomleft,~] = inpolygon(posrx,posry,roi_bottom_left.Position(:,1),roi_bottom_left.Position(:,2));
[bottomright,~] = inpolygon(posrx,posry,roi_bottom_right.Position(:,1),roi_bottom_right.Position(:,2));
[topright,~] = inpolygon(posrx,posry,roi_top_right.Position(:,1),roi_top_right.Position(:,2));

[topdoor,~] = inpolygon(posrx,posry,roi_top_door.Position(:,1),roi_top_door.Position(:,2));
[leftdoor,~] = inpolygon(posrx,posry,roi_left_door.Position(:,1),roi_left_door.Position(:,2));
[bottomdoor,~] = inpolygon(posrx,posry,roi_bottom_door.Position(:,1),roi_bottom_door.Position(:,2));
[rightdoor,~] = inpolygon(posrx,posry,roi_right_door.Position(:,1),roi_right_door.Position(:,2));



% figure;
% plot(posrx,posry,'color',[.6 .6 .6],'linewidth',2);
% hold on
% % plot position of position in each quadrant with different color
% scatter(posrx(topleft),posry(topleft),200,'r.'); % points inside
% scatter(posrx(bottomleft),posry(bottomleft),200,'g.'); % points inside
% scatter(posrx(bottomright),posry(bottomright),200,'b.'); 
% scatter(posrx(topright),posry(topright),200,'y.'); 
% scatter(posrx(topdoor),posry(topdoor),200,'.','MarkerEdgeColor',[1 0 1]); % points inside
% scatter(posrx(leftdoor),posry(leftdoor),200,'.','MarkerEdgeColor',[.3 .3 .3]); % points inside
% scatter(posrx(bottomdoor),posry(bottomdoor),200,'.','MarkerEdgeColor',[1 .5 0]); 
% scatter(posrx(rightdoor),posry(rightdoor),200,'.','MarkerEdgeColor',[.5 0 0]); 



