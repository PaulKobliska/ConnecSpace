Matlab code for ConnecSpace project.
This MATLAB code extract spatial and directionnal activity in individual rooms and entire in the entire environment.

Required Toolboxes : curve fitting, Image processing, signal processing & statistic and machine Learning 

Files
main_script.m: The main file to run the entire analysis.
functions folder: Contains helper functions used by the main script.
mat files: Folder containing example of a recording session data.

Usage
Open main_script.m in MATLAB.
Configure "pathfigure" and "results_folder" with paths where the figures and excel file will be saved.
Run the script
When a plot of the path of the animal appear it is needed to manualy delimit first, the left room 
(by clicking with the mouse and a double click to end a delimitation), then the right room and finally the door area.

Outputs
3 figures will be created:
1) the polar plot of the directionnal activity, the positions of the spikes superimposed with the trajectory of the animal, the circular autocorrelation and the bidirectional and tetradirectional scores.
2) the ratemaps of room 1, 2 and room 1 rotated to 180° and spatial correlations of unrotated rooms (1 vs 2) and rotated rooms (1 rotated vs 2).
3) Global ratemap, global polarplot, individual room ratemaps, individual rooms polarplots, door area ratemap and polarplot.

One excel file containing extracted datas.
