% This function allow to find the centre of the maze. It uses the
% coordinates of the maze measuered (they vary according to the recording
% room).
function centre = findCentre(NW,SW,NE,SE)

% The centre will be at the point of interception of the corner diagonals
a = (NE(2)-SW(2))/(NE(1)-SW(1)); % Slope for the NE-SW diagonal
b = (SE(2)-NW(2))/(SE(1)-NW(1)); % Slope for the SE-NW diagonal
c = SW(2);
d = NW(2);
x = (d-c+a*SW(1)-b*NW(1))/(a-b); % X-coord of centre
y = a*(x-SW(1))+c; % Y-coord of centre
centre = [x,y];