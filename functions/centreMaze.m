% This function allow to find the centre of the maze. It uses the
% coordinates of the maze measuered (they vary according to the recording
% room).
function [NE,NW,SW,SE] = centreMaze

% maze oriented as for rat SDL103!!!!!

% SET CORNERS OF THE MAZE (reference box that contains the maze)
%   NW = [80 230]; % value for Ingrid
%   SW = [80 30];
%   NE = [220 230];
%   SE = [220 30];

NW = [50 200]; % value for Sophie
SW = [50 20];
NE = [200 200];
SE = [200 20];