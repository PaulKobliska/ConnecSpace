% *************************************************************************
% script to analyze 2 and 4 compartments HD activities
% handmade script for 2019 M2 internship
% made by P-Yves 15/03/2019
%modification made for Celia PhD 05/01/2021
%modification made by PY 05/09/2022 : change selection of compartments and
%door by crating polygon with mouse clics.
% *************************************************************************

%define recording file (and tetrode and cell if wanted)
% inputfile = 'Germaine_12012021_S7-01.mat'; %recording file (.uff or .mat accepted)

% tetrode = 2; %[] if you want all cells orwrite the number if you want a specific cell
% cell = 1;   %[] write number of cell only if you specified a number of tetrode, else, write []
% compartment = 1 ; %set the number of compartment between 1, 2 or 4
% cd 'C:\Users\celia\Nextcloud\Shared\Celia script PhD\mat exemples';
% 
% cd 'C:\Users\celia\Pictures\ConnecSpace\FINAL ANALYSIS\Result\Protocole all cells\';
% [~,~,raw] = xlsread('Results_per_condition_Repeted_cells','S1 4 comp repeted');
% raw_number = str2double(raw);


cd 'C:\Users\celia\Documents\ConnecSpace\Analyses\Analyses ultimes\Mat\Protocole place cell\S1 4 iden - S2 4 diff - S3 4 diff+ 1 cue ext - S4 4 iden+ 1 cue ext\';
a= dir;
b= struct2table(a);
all_files= b.name(:);
% remove directories from files list
all_files([1,2],:)=[];
zz = size (all_files,1);

a=dir;  % Prendre tous les doc du fichier dans a
b=struct2table(a);  % convertir a en table
folder_name=b.name; % prendre seulement la colone "nom" des documents
folder_name(1:2)=[];    % enlever les 2 premières lignes de noms car elles sont vides

%define paramaters:
binWidth = 3;     % bining for ratemaps
smooth_factor = 2;  % factor of smoothing for ratemap
TimeWindow = 1000;  % length of time autocorrelogram (in msec)
Tbin = 10;          % bin of autocorrelogram (in msec)        
threshold = 100;    % for path correction
binSizeDir = 6;    % bining for polar plot (choose 5,6 or 10)
correction = 0;     % angles shift of LEDs position on the animal head
save_figure = 1;    % 1= save figure, 0=don't save figure (don't forget to set path where to save figure)
close_figure = 1;   % 1= close figure after display, 0=don't close figure
pathfigure = 'C:\Users\celia\Pictures\ConnecSpace\FINAL ANALYSIS\Images\'; % folder to save figure (finish the path with\)
results_folder = 'C:\Users\celia\Pictures\ConnecSpace\FINAL ANALYSIS\Result\';
mat_folder = 'C:\Users\celia\Pictures\ConnecSpace\FINAL ANALYSIS\Mat\';

%%
for f = 1%:zz

filename = all_files(f);
filename =filename{1};
inputfile = filename
% name = raw{f,1}; name = name(1:end-4);
% inputfile=strcat(name,'_t',num2str(raw{f,2}),'_c',num2str(raw{f,3}),'.mat');

% *************************************************************************
shape = 'circle';

if strcmp(inputfile(end-3:end),'.uff')
[EEGvalues, EEGts, TTMtx, PosMtx, E_FlagMtx, eeg_chan_to_read,inputfile] = Read_UFF_Slow(inputfile, 1, 1, 0, 0, 0);

elseif strcmp(inputfile(end-3:end),'.mat')
load(inputfile);
end

[posx,posy,posts,mapAxis,visited] = posdata_multicompartment(PosMtx,threshold,shape,binWidth);
path = figure; plot(posx,posy,'color',[.6 .6 .6]); %hold on; plot(posgx,posgy,'color',[.2 .2 .2]);

answer = questdlg('How many compartment(s)?','Number of compartment(s)', ...
    '1','2','4','4');
% Handle response
switch answer
    case '1'
        compartment =1;
    case '2'
        compartment =2;
    case '4'
        compartment =4;
end
close(path);

    if compartment == 2
%         if strcmp(raw(f,3),raw(f-1,3))==0 %folder_name{m}(1:end-4),folder_name{m-1}(1:end-4)))==0
%         f=[];raw=[];
        [xy_left,xy_right,xy_door,left,right,door] = Splitting_compartment_v1_removing_door(PosMtx,inputfile);
        close
%         end
    end
    
        if compartment == 4
%             if strcmp(raw(f,3),raw(f-1,3))==0
               [xy_top_left,xy_bottom_left,xy_top_right,xy_bottom_right,xy_top_door,xy_left_door,xy_bottom_door,xy_right_door,topleft,bottomleft,bottomright,topright,topdoor,leftdoor,bottomdoor,rightdoor] = Splitting_four_compartment_pyves_removing_door(PosMtx);
               close
%               end
        end
        TTMtx=TTMtx2;
 
tetrode = unique(TTMtx(:,2));    
    for pp = 2%:length(tetrode)
    TTMtx2=TTMtx(TTMtx(:,2)== tetrode(pp),:);
    TTMtx2(TTMtx2(:,3)==0,:) = [];
    if isempty(TTMtx2)
    else 
   p = tetrode(pp);

    cell = unique(TTMtx2(:,3));
    for iii = 1:length(cell)
    ts = TTMtx2(TTMtx2(:,3)==cell(iii),1);
    cells = cell(iii);
    waveforms = TTMtx2(TTMtx2(:,3)==cell(iii),4:131);
    i = cells;

    %% check whether waveform if correct
%    [selected_el,average_selected_el,std_plus_selected_el,std__selected_el,std_minus_selected_el,...
%     average_derive_selected_el,std__derive_selected_el,std_plus_derive_selected_el,std_minus_derive_selected_el]= waveform(waveforms);
%    waveformtocheck = figure;
%    plot(average_selected_el,'linewidth',3,'color','k'); axis square;
%     hold on
%     plot(std_plus_selected_el,'--','linewidth',3,'color','k'); axis square;
%     plot(std_minus_selected_el,'--','linewidth',3,'color','k'); axis square;
%     hold off
% %     axis off
%     box off
%     title('Mean waveform','fontsize',12);
%     xlim([0 32]);
%     set(gca, 'XTickLabel', [])
%     set(gca, 'XTick', [])
%     set(gca, 'YTickLabel', [])
%     set(gca, 'YTick', [])
%     h=get(gcf, 'currentaxes');
%     set(h, 'fontsize', 12, 'linewidth', 2);

%     answer = questdlg('looks like','Waveform?', ...
%     'a neuron','a bullshit','a bullshit');
% %     handle response
% switch answer
%     case 'a neuron'
%         visualcheckingofwaveform =1;
%     case 'a bullshit'
%         visualcheckingofwaveform =0;
% end
% close(waveformtocheck);
%    
%     if visualcheckingofwaveform == 1

if compartment == 1
%%
general_fig_1comp = figure('visible', 'off');
[selected_el,average_selected_el,std_plus_selected_el,std__selected_el,std_minus_selected_el,...
    average_derive_selected_el,std__derive_selected_el,std_plus_derive_selected_el,std_minus_derive_selected_el]= waveform(waveforms);
subplot(2,3,1);
    plot(average_selected_el,'linewidth',3,'color','k'); axis square;
    hold on
    plot(std_plus_selected_el,'--','linewidth',3,'color','k'); axis square;
    plot(std_minus_selected_el,'--','linewidth',3,'color','k'); axis square;
    hold off
    axis off
    box off
    title('Mean waveform','fontsize',12);
    xlim([0 32]);
    set(gca, 'XTickLabel', [])
    set(gca, 'XTick', [])
    set(gca, 'YTickLabel', [])
    set(gca, 'YTick', [])
    h=get(gcf, 'currentaxes');
    set(h, 'fontsize', 12, 'linewidth', 2);

% subplot(2,3,2);
% plot(posx,posy,'color',[.4 .4 .4],'linewidth',2);%,spkx,spky,'.r');
% hold on
% [spkx,spky] = get_pos_spikes(ts,posx,posy,posts); % get the position of the spikes
% scatter(spkx,spky,400,'.r');
% set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);
% axis off
% box off

subplot(2,3,2);
plot(posx,posy,'color',[.4 .4 .4],'linewidth',2);
% patch(posrx,posry,[.2 .2 .2],'EdgeColor',[.2 .2 .2],'LineWidth',2,'FaceColor','non');
[mapAxis] = mapaxis(posx,posy,binWidth); visited = visitedBins(posx,posy,mapAxis);
[spkx,spky] = get_pos_spikes(ts,posx,posy,posts); % get the position of the spikes
hold on
scatter(spkx,spky,400,'.r');
set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);
axis off

[ratemap] = ratemap_gaussian(2*binWidth,spkx,spky,posx,posy,posts,binWidth,mapAxis); % build matrix ratemap 
ratemap(visited==0) = NaN;
subplot(2,3,3);
myfig=pcolor(ratemap);
colormap jet
set(myfig,'Edgecolor','none')
axis xy
set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);
axis off
box off
peak_rate=max(ratemap(:))*1000; peak_rate_comp=peak_rate;peak_rate_comp=sprintf('%.2f', peak_rate_comp);
title(strcat('Peak rate = ',peak_rate_comp, ' Hz'));
peak_rate_global = peak_rate;

% subplot(2,3,3);
% 
% polarplot_compartment_new(posts,posgx,posgy,posrx,posry,ts,correction,binSizeDir);
% set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);

% calculate and plot cell firing/speed relationship
subplot(2,3,4);
[ret,beta,f0]=speed_firing(posx,posy,posts,ts);
h=get(gcf, 'currentaxes');
set(h, 'fontsize', 15, 'linewidth', 2);
box off;
get(gca,'ticklength');
set(gca,'TickLength',[0.03,0.02]);
set(gca,'TickDir','out');
    
% calculate and plot time firing autocorrelagram
subplot(2,3,5);
[autocorr,y] = autocorrelogram(TimeWindow,Tbin,ts);
h=get(gcf, 'currentaxes');
set(h, 'fontsize', 15, 'linewidth', 2);
box off;
get(gca,'ticklength');
set(gca,'TickLength',[0.03,0.02]);
set(gca,'TickDir','out');
title('Autocorrelogram','fontsize',12);

subplot(2,3,6);
[mean_freq,freq_peak_distr,peak_theta,peak_delta] = intrinsicfiring(TimeWindow,Tbin,autocorr,y);
title(strcat('Peak theta=',num2str(sprintf('%.2f', freq_peak_distr),'Hz')),'fontsize',12);
[peak_frequency_theta,theta_power,delta_power,sig_threshold,theta_modulation] = theta_delta(TimeWindow,Tbin,ts);
set(gca,'TickDir','out','Box','off');
h=get(gcf, 'currentaxes');
set(h, 'fontsize', 12, 'linewidth', 2);
set(gca,'TickLength',[0.02,0.01]);
title(strcat('Peak theta=',num2str(sprintf('%.2f', freq_peak_distr),'Hz')),'fontsize',12);


%% save general figure
set(general_fig_1comp,'units','normalized','outerposition',[0 0 1 1 ]);
set(general_fig_1comp,'Units','Inches');
pos = get(general_fig_1comp,'Position');
set(general_fig_1comp,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

   if save_figure == 1
            complete_name = strcat(pathfigure,inputfile(1:end-4),'_t',num2str(p),'_c',num2str(i),'_general_one_comp');
            print('-dpng',complete_name,'-r100');
%             print(two_comp,complete_name,'-dpdf')
   end
    if close_figure == 1
        close(general_fig_1comp);
    end

%% colored polar plot and trajectory X spikes map
posts = PosMtx(:,1);
posgx = PosMtx(:,2);
posgy = PosMtx(:,3);
[posgx,posgy] = smooth_path(posgx,posgy);
posts = PosMtx(:,1);
posrx = PosMtx(:,4);
posry = PosMtx(:,5);
[posrx,posry] = smooth_path(posrx,posry);
zztop = figure('visible', 'off');
[PFD_deg,Vector_length,peak_HD,~,deg_theta,posdirts,dirmap]= polarplot_compartment_new(posts,posgx,posgy,posrx,posry,ts,correction,binSizeDir);
close(zztop);

colored_plot = figure('visible', 'off');
[dirrate] = polarplot_v2_colored_new(PosMtx,spkx,spky,ts,binSizeDir,correction,shape);

set(colored_plot,'units','normalized','outerposition',[0 0 1 1]);
set(colored_plot,'Units','Inches');
pos = get(colored_plot,'Position');
set(colored_plot,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

[Vector_length_double,flipping_score,flipping_score_range,Vector_length_quadri,Four_dir_score,Four_dir_score_minmax,Four_dir_score_range,Four_dir_score_minmax_range,correlation,minimums,maximums,corr_value,corr_bin]=flipping_scores(dirrate,binSizeDir,PFD_deg);
minimums_angles = corr_bin(minimums);
maximums_angles = corr_bin(maximums);
minimums_values = corr_value(minimums);
maximums_values = corr_value(maximums);

number_of_corr_peaks = sum(minimums)+sum(maximums);

subplot(2,2,3);
plot(correlation(:,2),correlation(:,1),'k','linewidth',2); hold on;
scatter(corr_bin(maximums),corr_value(maximums),100,'o','markeredgecolor',[.3 .3 .3],'markerfacecolor',[.7 .7 .7]);
scatter(corr_bin(minimums),corr_value(minimums),100,'o','markeredgecolor',[.3 .3 .3],'markerfacecolor',[.7 .7 .7]);
axis square; box off;
h=get(gcf, 'currentaxes');
set(h, 'fontsize', 15, 'linewidth', 2);
box off;
get(gca,'ticklength');
set(gca,'TickLength',[0.03,0.02]);
set(gca,'TickDir','out');
xlim([0 360]);
xticks(0:45:360);
xlabel('degre'); ylabel('correlation coefficient');

subplot(2,2,4); axis off;
text(0,1,strcat('flip score = ',num2str(flipping_score)), 'fontsize', 20);
text(0,0.8,strcat('flips score range = ',num2str(flipping_score_range)), 'fontsize', 20);
text(0,0.6,strcat('quadiscore = ',num2str(Four_dir_score)), 'fontsize', 20);
text(0,0.4,strcat('quadiscore range = ',num2str(Four_dir_score_range)), 'fontsize', 20);
text(0,0.2,strcat('quadiscore minmax = ',num2str(Four_dir_score_minmax)), 'fontsize', 20);
text(0,0,strcat('quadiscore minmax range = ',num2str(Four_dir_score_minmax_range)), 'fontsize', 20);

   if save_figure == 1
            complete_name = strcat(pathfigure,inputfile(1:end-4),'_t',num2str(p),'_c',num2str(i),'_colored_plot_two_comp');
            print('-dpng',complete_name,'-r500');
%             print(colored_plot,complete_name,'-dpdf');

   end
   if close_figure == 1
        close(colored_plot);
   end
   
%% save results
results_txt_1comp(results_folder,inputfile,p,i,peak_rate_global,ret,...
    PFD_deg,Vector_length,peak_HD,...
    Vector_length_double,flipping_score,flipping_score_range,...
    Vector_length_quadri,Four_dir_score,Four_dir_score_minmax,Four_dir_score_range,Four_dir_score_minmax_range,...
    peak_frequency_theta,theta_power,delta_power,sig_threshold,theta_modulation,number_of_corr_peaks)

%% save mat
nameformat = strcat(mat_folder,filename(1:end-4),'_t',num2str(p),'_c',num2str(i),'.mat');
save(nameformat,...
    'TTMtx2','PosMtx','filename','ts','waveforms','ratemap','dirrate');
        
%% Splitting 2 compartments
elseif compartment == 2
   [spkx,spky,spike_left,spike_right,spike_door] = Splitting_two_compartment_and_door(ts,posx,posy,posts,xy_left,xy_right,xy_door);
%%
two_comp = figure('visible', 'off');
 
%% Rate map global
% % % [posx,posy] = smooth_path(posx,posy); [posx,posy] = center_path(posx,posy,shape);
% % % [mapAxis] = mapaxis(posx,posy,binWidth); visited = visitedBins(posx,posy,mapAxis);   
[spkx,spky] = get_pos_spikes(ts,posx,posy,posts); % get the position of the spikes
[ratemap] = ratemap_gaussian(2*binWidth,spkx,spky,posx,posy,posts,binWidth,mapAxis); % build matrix ratemap 
ratemap(visited==0) = NaN;
subplot(4,3,1);
myfig=pcolor(ratemap);
colormap jet
set(myfig,'Edgecolor','none')
axis xy
set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);
axis off
box off
peak_rate=max(ratemap(:))*1000; peak_rate_global=peak_rate;peak_rate_global=sprintf('%.2f', peak_rate_global);
title(strcat('Peak rate = ',peak_rate_global, ' Hz'));
peak_rate_global = peak_rate;

%% polar plot global
posts = PosMtx(:,1);
posgx = PosMtx(:,2);
posgy = PosMtx(:,3);
[posgx,posgy] = smooth_path(posgx,posgy);
posts = PosMtx(:,1);
posrx = PosMtx(:,4);
posry = PosMtx(:,5);
[posrx,posry] = smooth_path(posrx,posry);

subplot(4,3,2);
[PFD_deg,Vector_length,peak_HD,dirrate,deg_theta,posdirts,dirmap]= polarplot_compartment_new(posts,posgx,posgy,posrx,posry,ts,correction,binSizeDir);
peak_HD_2=sprintf('%.2f', peak_HD);
Vector_length_2=sprintf('%.2f', Vector_length);
set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);
title({
    ['Pic = ', peak_HD_2,'Hz, Vector length = ',Vector_length_2]},'fontsize',12);
set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);xlabel('global')


%% polar plots for 2 compartments and calculate PFD, vector length and peak rate   

%path for each compartment and the door
subplot(4,3,3);
plot(posrx,posry);
% patch(posrx,posry,[.2 .2 .2],'EdgeColor',[.2 .2 .2],'LineWidth',2,'FaceColor','non');

hold on
% plot position of position in each quadrant with different color
scatter(posrx(left),posry(left),'r.'); % points inside
scatter(posrx(right),posry(right),'g.'); % points inside
scatter(posrx(door),posry(door),'b.'); 
axis square
axis off

subplot(4,3,10);
% polar plot left compartment
[PFD_deg_left,Vector_length_left,peak_HD_left,dirrate_left,deg_theta,posdirts,dirmap]= polarplot_compartment_new(posts(left),posgx(left),posgy(left),posrx(left),posry(left),ts(spike_left),correction,binSizeDir);
peak_HD_2=sprintf('%.2f', peak_HD_left);
Vector_length_2=sprintf('%.2f', Vector_length_left);
set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);
title({['Pic = ',peak_HD_2,'Hz, Vector length=',Vector_length_2]});
set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);

subplot(4,3,11);
% polar plot right compartment
[PFD_deg_right,Vector_length_right,peak_HD_right,dirrate_right,deg_theta,posdirts,dirmap]= polarplot_compartment_new(posts(right),posgx(right),posgy(right),posrx(right),posry(right),ts(spike_right),correction,binSizeDir);
peak_HD_2=sprintf('%.2f', peak_HD_right);
Vector_length_2=sprintf('%.2f', Vector_length_right);
set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);
title({['Pic = ',peak_HD_2,'Hz, Vector length=',Vector_length_2]});
set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);

subplot(4,3,12);
% polar plot door
[PFD_deg_door,Vector_length_door,peak_HD_door,dirrate_door,deg_theta,posdirts,dirmap]= polarplot_compartment_new(posts(door),posgx(door),posgy(door),posrx(door),posry(door),ts(spike_door),correction,binSizeDir);
peak_HD_2=sprintf('%.2f', peak_HD_door);
Vector_length_2=sprintf('%.2f', Vector_length_door);
set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);
title({['Pic = ',peak_HD_2,'Hz, Vector length=',Vector_length_2]});
set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);

%%
% ratemap left compartment
[posrx,posry] = center_path(posrx,posry,shape);

subplot(4,3,4);
plot(posrx,posry,'color',[.2 .2 .2],'linewidth',2);
% patch(posrx,posry,[.2 .2 .2],'EdgeColor',[.2 .2 .2],'LineWidth',2,'FaceColor','non');
[mapAxis] = mapaxis(posrx,posry,binWidth); visited = visitedBins(posrx(left),posry(left),mapAxis);
ts_left = ts(spike_left);
[spkx_left,spky_left] = get_pos_spikes(ts_left,posrx(left),posry(left),posts(left)); % get the position of the spikes
hold on
scatter(spkx_left,spky_left,200,'.r');
set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);
axis off

[ratemap_left] = ratemap_gaussian(2*binWidth,spkx_left,spky_left,posrx(left),posry(left),posts(left),binWidth,mapAxis); % build matrix ratemap 
ratemap_left(visited==0) = NaN;
subplot(4,3,7);
myfig=pcolor(ratemap_left);
colormap jet
set(myfig,'Edgecolor','none')
axis xy
set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);
axis off
box off
peak_rate=max(ratemap_left(:))*1000; peak_rate_left=peak_rate;peak_rate_left=sprintf('%.2f', peak_rate_left);
title(strcat('Peak rate = ',peak_rate_left, ' Hz'));
peak_rate_left = peak_rate;

% ratemap right compartment

subplot(4,3,5);
plot(posrx,posry,'color',[.2 .2 .2],'linewidth',2);
% patch(posrx,posry,[.2 .2 .2],'EdgeColor',[.2 .2 .2],'LineWidth',2,'FaceColor','non');

[mapAxis] = mapaxis(posrx,posry,binWidth); visited = visitedBins(posrx(right),posry(right),mapAxis);
ts_right = ts(spike_right);
[spkx_right,spky_right] = get_pos_spikes(ts_right,posrx(right),posry(right),posts(right)); % get the position of the spikes
hold on
scatter(spkx_right,spky_right,200,'.r');
set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);
axis off

[ratemap_right] = ratemap_gaussian(2*binWidth,spkx_right,spky_right,posrx(right),posry(right),posts(right),binWidth,mapAxis); % build matrix ratemap 
ratemap_right(visited==0) = NaN;
subplot(4,3,8);
myfig=pcolor(ratemap_right);
colormap jet
set(myfig,'Edgecolor','none')
axis xy
set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);
axis off
box off
peak_rate=max(ratemap_right(:))*1000; peak_rate_right=peak_rate;peak_rate_right=sprintf('%.2f', peak_rate_right);
title(strcat('Peak rate = ',peak_rate_right, ' Hz'));
peak_rate_right = peak_rate;

subplot(4,3,6); % Door
plot(posrx,posry,'color',[.2 .2 .2],'linewidth',2);
% patch(posrx,posry,[.2 .2 .2],'EdgeColor',[.2 .2 .2],'LineWidth',2,'FaceColor','non');

[mapAxis] = mapaxis(posrx,posry,binWidth); visited = visitedBins(posrx,posry,mapAxis);
ts_door = ts(spike_door);
[spkx_door,spky_door] = get_pos_spikes(ts_door,posrx(door),posry(door),posts(door)); % get the position of the spikes
hold on
scatter(spkx_door,spky_door,200,'.r');
set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);
axis off
hold off

[ratemap_door] = ratemap_gaussian(2*binWidth,spkx_door,spky_door,posrx(door),posry(door),posts(door),binWidth,mapAxis); % build matrix ratemap 
ratemap_door(visited==0) = NaN;
subplot(4,3,9);
myfig=pcolor(ratemap_door);
colormap jet
set(myfig,'Edgecolor','none')
axis xy
set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);
axis off
box off
peak_rate=max(ratemap_door(:))*1000; peak_rate_door=peak_rate;peak_rate_door=sprintf('%.2f', peak_rate_door);
title(strcat('Peak rate = ',peak_rate_door, ' Hz'));
peak_rate_door = peak_rate;

%% save figures of 2 compartments polar plots  
set(two_comp,'units','normalized','outerposition',[0 0 1 1 ]);
set(two_comp,'Units','Inches');
pos = get(two_comp,'Position');
set(two_comp,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

   if save_figure == 1
            complete_name = strcat(pathfigure,inputfile(1:end-4),'_t',num2str(p),'_c',num2str(i),'_split_two_comp');
            print('-dpng',complete_name,'-r100');
%             print(two_comp,complete_name,'-dpdf')
   end
    if close_figure == 1
        close(two_comp);
    end


%%
general_fig = figure('visible', 'off');
% [selected_el,average_selected_el,std_plus_selected_el,std__selected_el,std_minus_selected_el,...
%     average_derive_selected_el,std__derive_selected_el,std_plus_derive_selected_el,std_minus_derive_selected_el]= waveform(waveforms);
subplot(2,3,1);
    plot(average_selected_el,'linewidth',3,'color','k'); axis square;
    hold on
    plot(std_plus_selected_el,'--','linewidth',3,'color','k'); axis square;
    plot(std_minus_selected_el,'--','linewidth',3,'color','k'); axis square;
    hold off
    axis off
    box off
    title('Mean waveform','fontsize',12);
    xlim([0 32]);
    set(gca, 'XTickLabel', [])
    set(gca, 'XTick', [])
    set(gca, 'YTickLabel', [])
    set(gca, 'YTick', [])
    h=get(gcf, 'currentaxes');
    set(h, 'fontsize', 12, 'linewidth', 2);

subplot(2,3,2);
plot(posx,posy,'color',[.4 .4 .4],'linewidth',2);%,spkx,spky,'.r');
 hold on
scatter(spkx,spky,400,'.r');
set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);
axis off
box off

posts = PosMtx(:,1);
posgx = PosMtx(:,2);
posgy = PosMtx(:,3);
[posgx,posgy] = smooth_path(posgx,posgy);
posts = PosMtx(:,1);
posrx = PosMtx(:,4);
posry = PosMtx(:,5);
[posrx,posry] = smooth_path(posrx,posry);
subplot(2,3,3);
polarplot_compartment_new(posts,posgx,posgy,posrx,posry,ts,correction,binSizeDir);
set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);

% calculate and plot cell firing/speed relationship
subplot(2,3,4);
[ret,beta,f0]=speed_firing(posx,posy,posts,ts);
h=get(gcf, 'currentaxes');
set(h, 'fontsize', 15, 'linewidth', 2);
box off;
get(gca,'ticklength');
set(gca,'TickLength',[0.03,0.02]);
set(gca,'TickDir','out');
    
% calculate and plot time firing autocorrelagram
subplot(2,3,5);
[autocorr,y] = autocorrelogram(TimeWindow,Tbin,ts);
h=get(gcf, 'currentaxes');
set(h, 'fontsize', 15, 'linewidth', 2);
box off;
get(gca,'ticklength');
set(gca,'TickLength',[0.03,0.02]);
set(gca,'TickDir','out');
title('Autocorrelogram','fontsize',12);

subplot(2,3,6);
[mean_freq,freq_peak_distr,peak_theta,peak_delta] = intrinsicfiring(TimeWindow,Tbin,autocorr,y);
title(strcat('Peak theta=',num2str(sprintf('%.2f', freq_peak_distr),'Hz')),'fontsize',12);
[peak_frequency_theta,theta_power,delta_power,sig_threshold,theta_modulation] = theta_delta(TimeWindow,Tbin,ts);
set(gca,'TickDir','out','Box','off');
h=get(gcf, 'currentaxes');
set(h, 'fontsize', 12, 'linewidth', 2);
set(gca,'TickLength',[0.02,0.01]);
title(strcat('Peak theta=',num2str(sprintf('%.2f', freq_peak_distr),'Hz')),'fontsize',12);


%% save general figures 
set(general_fig,'units','normalized','outerposition',[0 0 1 1 ]);
set(general_fig,'Units','Inches');
pos = get(general_fig,'Position');
set(general_fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

   if save_figure == 1
            complete_name = strcat(pathfigure,inputfile(1:end-4),'_t',num2str(p),'_c',num2str(i),'_general_two_comp');
            print('-dpng',complete_name,'-r100');
%             print(two_comp,complete_name,'-dpdf')
   end
    if close_figure == 1
        close(general_fig);
    end

%% ratemap comp 1
[posrx_comp1,posry_comp1] = center_path(posrx(left),posry(left),shape);
[mapAxis_comp1] = mapaxis(posrx_comp1,posry_comp1,binWidth);
visited_comp1 = visitedBins(posrx_comp1,posry_comp1,mapAxis_comp1); % calculate the visited zones of the arena 
[spkx_comp1,spky_comp1] = get_pos_spikes(ts(spike_left),posrx_comp1,posry_comp1,posts(left)); % get the position of the spikes
[ratemap_comp1] = ratemap_gaussian(2*binWidth,spkx_comp1,spky_comp1,posrx_comp1,posry_comp1,posts(left),binWidth,mapAxis_comp1); % build matrix ratemap 
ratemap_comp1(visited_comp1==0) = NaN;

% ratemap comp 1 rotated by 180°
minus_posrx_comp1 = -(posrx_comp1);
minus_posry_comp1 = -(posry_comp1);
[minus_posrx_comp1,minus_posry_comp1] = center_path(minus_posrx_comp1,minus_posry_comp1,shape);
[rotated_mapAxis_comp1] = mapAxis_comp1;%mapaxis(minus_posrx_comp1,minus_posry_comp1,binWidth);
rotated_visited_comp1 = visitedBins(minus_posrx_comp1,minus_posry_comp1,rotated_mapAxis_comp1); % calculate the visited zones of the arena 
[rotated_spkx_comp1,rotated_spky_comp1] = get_pos_spikes(ts(spike_left),minus_posrx_comp1,minus_posry_comp1,posts(left)); % get the position of the spikes
[rotated_ratemap_comp1] = ratemap_gaussian(2*binWidth,rotated_spkx_comp1,rotated_spky_comp1,minus_posrx_comp1,minus_posry_comp1,posts(left),binWidth,rotated_mapAxis_comp1); % build matrix ratemap 
rotated_ratemap_comp1(rotated_visited_comp1==0) = NaN;

[posrx_comp2,posry_comp2] = center_path(posrx(right),posry(right),shape);

[mapAxis_comp2] = mapAxis_comp1;% mapaxis(posrx_comp2,posry_comp2,binWidth);
visited_comp2 = visitedBins(posrx_comp2,posry_comp2,mapAxis_comp2); % calculate the visited zones of the arena 
[spkx_comp2,spky_comp2] = get_pos_spikes(ts(spike_right),posrx_comp2,posry_comp2,posts(right)); % get the position of the spikes

[ratemap_comp2] = ratemap_gaussian(2*binWidth,spkx_comp2,spky_comp2,posrx_comp2,posry_comp2,posts(right),binWidth,mapAxis_comp2); % build matrix ratemap 
    % remove invisited bins from the rate map
ratemap_comp2(visited_comp2==0) = NaN; 

%% Compartment correlation
[corrcoeff,nonZero_corrcoeff,RrankO] = get_2D_correlation_value(ratemap_comp1,ratemap_comp2);
[rotated_corrcoeff,rotated_nonZero_corrcoeff,rotated_RrankO] = get_2D_correlation_value(rotated_ratemap_comp1,ratemap_comp2);

%% save figures of 2 compartments corr  

two_comp_corr = figure('visible', 'off');
subplot(2,2,1);
myfig=pcolor(ratemap_comp1);
colormap jet
set(myfig,'Edgecolor','none')
axis xy
set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);
axis off
box off
title('compartiment 1');
subplot(2,2,2);
myfig=pcolor(ratemap_comp2);
colormap jet
set(myfig,'Edgecolor','none')
axis xy
set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);
axis off
box off
title('compartiment 2');

subplot(2,2,3);
myfig=pcolor(rotated_ratemap_comp1);
colormap jet
set(myfig,'Edgecolor','none')
axis xy
set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);
axis off
box off
title('rotated compartiment 1');
subplot(2,2,4);
text(0,1,strcat('corr 0° =',num2str (corrcoeff)),'FontSize', 30);
text(0,0.5,strcat('corr 180° =',num2str (rotated_corrcoeff)),'FontSize', 30);
axis off

set(two_comp_corr,'units','normalized','outerposition',[0 0 1 1 ]);
set(two_comp_corr,'Units','Inches');
pos = get(two_comp_corr,'Position');
set(two_comp_corr,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

   if save_figure == 1
            complete_name = strcat(pathfigure,inputfile(1:end-4),'_t',num2str(p),'_c',num2str(i),'_spatial_correlation_two_comp');
            print('-dpng',complete_name,'-r100');
%             print(two_comp_corr,complete_name,'-dpdf')

   end
   if close_figure == 1
        close(two_comp_corr);
   end

%% colored polar plot and trajectory X spikes map

colored_plot = figure('visible', 'off');
[dirrate] = polarplot_v2_colored_new(PosMtx,spkx,spky,ts,binSizeDir,correction,shape);

set(colored_plot,'units','normalized','outerposition',[0 0 1 1 ]);
set(colored_plot,'Units','Inches');
pos = get(colored_plot,'Position');
set(colored_plot,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

[Vector_length_double,flipping_score,flipping_score_range,Vector_length_quadri,Four_dir_score,Four_dir_score_minmax,Four_dir_score_range,Four_dir_score_minmax_range,correlation,minimums,maximums,corr_value,corr_bin]=flipping_scores(dirrate,binSizeDir,PFD_deg);
minimums_angles = corr_bin(minimums);
maximums_angles = corr_bin(maximums);
minimums_values = corr_value(minimums);
maximums_values = corr_value(maximums);

number_of_corr_peaks = sum(minimums)+sum(maximums);


subplot(2,2,3);
plot(correlation(:,2),correlation(:,1),'k','linewidth',2); hold on;
scatter(corr_bin(maximums),corr_value(maximums),100,'o','markeredgecolor',[.3 .3 .3],'markerfacecolor',[.7 .7 .7]);
scatter(corr_bin(minimums),corr_value(minimums),100,'o','markeredgecolor',[.3 .3 .3],'markerfacecolor',[.7 .7 .7]);
axis square; box off;
h=get(gcf, 'currentaxes');
set(h, 'fontsize', 15, 'linewidth', 2);
box off;
get(gca,'ticklength');
set(gca,'TickLength',[0.03,0.02]);
set(gca,'TickDir','out');
xlim([0 360]);
xticks(0:45:360);
xlabel('degre'); ylabel('correlation coefficient');

subplot(2,2,4); axis off;
text(0,1,strcat('flip score = ',num2str(flipping_score)), 'fontsize', 20);
text(0,0.8,strcat('flips score range = ',num2str(flipping_score_range)), 'fontsize', 20);
text(0,0.6,strcat('quadiscore = ',num2str(Four_dir_score)), 'fontsize', 20);
text(0,0.4,strcat('quadiscore range = ',num2str(Four_dir_score_range)), 'fontsize', 20);
text(0,0.2,strcat('quadiscore minmax = ',num2str(Four_dir_score_minmax)), 'fontsize', 20);
text(0,0,strcat('quadiscore minmax range = ',num2str(Four_dir_score_minmax_range)), 'fontsize', 20);


   if save_figure == 1
            complete_name = strcat(pathfigure,inputfile(1:end-4),'_t',num2str(p),'_c',num2str(i),'_colored_plot_two_comp');
            print('-dpng',complete_name,'-r100');
%             print(colored_plot,complete_name,'-dpdf');

   end
   if close_figure == 1
        close(colored_plot);
   end
%% temporal analysis of cell activity per compartment
posts_left = posts(left);

diffposts_in = diff(posts_left);
diffsuperioto500ms = find(diffposts_in>500);


rate_in = [];duration_in = [];
for hgt = 1:length(diffsuperioto500ms)-1
%     path_x = posrx_in(diffsuperioto500ms(hgt)+1:diffsuperioto500ms(hgt+1));
%     path_y = posry_in(diffsuperioto500ms(hgt)+1:diffsuperioto500ms(hgt+1));
    time = posts_left(diffsuperioto500ms(hgt)+1:diffsuperioto500ms(hgt+1));
    spike = find(ts(spike_left)>posts_left(diffsuperioto500ms(hgt)+1) & ts(spike_left)<posts_left(diffsuperioto500ms(hgt+1)));
    duration = (time(end)-time(1))/1000;
    numspike = length(spike);
    rate_in(hgt)=numspike/duration;
    duration_in(hgt) = duration;
    
end

rate_in_without0= rate_in;
rate_in_without0(rate_in_without0==0) = [];


% rate_in(duration_in<1)=[];
% duration_in(duration_in<1)=[];
% duration_in(rate_in==0) = [];
% rate_in(rate_in==0) = [];

posts_right = posts(right);
diffposts_out = diff(posts_right);
   diffsuperioto500ms = find(diffposts_out>500);

rate_out = [];duration_out = [];
for hgt = 1:length(diffsuperioto500ms)-1
    time = posts_right(diffsuperioto500ms(hgt)+1:diffsuperioto500ms(hgt+1));
    spike = find(ts(spike_right)>posts_right(diffsuperioto500ms(hgt)+1) & ts(spike_right)<posts_right(diffsuperioto500ms(hgt+1)));
    duration = (time(end)-time(1))/1000;
    numspike = length(spike);
    rate_out(hgt)=numspike/duration;
    duration_out(hgt) = duration;
end
% rate_out(duration_out<1)=[];
% duration_out(duration_out<1)=[];
% duration_out(rate_out==0) = [];
% rate_out(rate_out==0) = [];
rate_out_without0= rate_out;
rate_out_without0(rate_out_without0==0) = [];


%% raw temporal analysis of cell activity (compartment colored)
path_order = zeros(length(posts),1);
path_order(left)=10;
path_order(right)=-10;
path_order(door)=1;
path_order(path_order==0)=1000;

if path_order(end)==1000
        gh = length(path_order);
        while path_order(gh) == 1000
            path_order(gh) = [];
            gh = gh-1;
        end    
end

copy_path = path_order;
for ii = 1:length(path_order)

    if path_order(ii) == 1000
        j = ii;
        while path_order(j) == 1000
            j = j+1;
        end

        next_val = path_order(j);
        copy_path(ii) = next_val;


    end
end

transition = abs(diff(copy_path));
[id,~] = find(transition>0);

remove_short_id = diff(id);
id(remove_short_id==1)=[];

rate = [];duration_path = [];

    time = posts(1:id(1));
    spike = find(ts>posts(1) & ts<posts(id(1)));
    duration = (time(end)-time(1))/1000;
    numspike = length(spike);
    rate(1,1)=numspike/duration;
    duration_path(1,1) = duration;

for thtg = 1:length(id)-1
    time = posts(id(thtg):id(thtg+1));
    spike = find(ts>posts(id(thtg)) & ts<posts(id(thtg+1)));
    duration = (time(end)-time(1))/1000;
    numspike = length(spike);
    rate(thtg+1)=numspike/duration;
    duration_path(thtg+1) = duration;   
end

    time = posts(id(end):end);
    spike = find(ts>posts(id(end)) & ts<posts(end));
    duration = (time(end)-time(1))/1000;
    numspike = length(spike);
    rate(1,end+1)=numspike/duration;
    duration_path(1,end+1) = duration;

compartment_side = copy_path(id);
rate = rate';
rate(rate == 0) = -1;
idx_left = find(compartment_side == 10);
idx_right = find(compartment_side == -10);
idx_door = find(compartment_side == 1);
  
%% plots temporal analysis for 2 compartments  
temporal_fig = figure('visible', 'off');

[posrx,posry] = center_path(posrx,posry,shape);

subplot(3,3,1);
plot(posrx,posry,'color',[.2 .2 .2],'linewidth',2);
 hold on
scatter(spkx_left,spky_left,200,'.r');
set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);
axis off
box off

subplot(3,3,4);
plot(posrx,posry,'color',[.2 .2 .2],'linewidth',2);
 hold on
scatter(spkx_right,spky_right,200,'.r');
set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);
axis off
box off

subplot(3,3,2);
plot(rate_in,'k-o','LineWidth',1,'MarkerSize',5,'MarkerEdgeColor','k','MarkerFaceColor',[0.5,0.5,0.5]);
axis square; box off;
h=get(gcf, 'currentaxes');
set(h, 'fontsize', 15, 'linewidth', 2);
box off;
get(gca,'ticklength');
set(gca,'TickLength',[0.03,0.02]);
set(gca,'TickDir','out');
ylim([0 max([rate_in,rate_out])]);
xlabel('path'); ylabel('frequency(Hz)');

subplot(3,3,5);
plot(rate_out,'k-o','LineWidth',1,'MarkerSize',5,'MarkerEdgeColor','k','MarkerFaceColor',[0.5,0.5,0.5]);
axis square; box off;
h=get(gcf, 'currentaxes');
set(h, 'fontsize', 15, 'linewidth', 2);
box off;
get(gca,'ticklength');
set(gca,'TickLength',[0.03,0.02]);
set(gca,'TickDir','out');
ylim([0 max([rate_in,rate_out])]);
xlabel('path'); ylabel('frequency(Hz)');

subplot(3,3,3);
plot(rate_in_without0,'k-o','LineWidth',2,'MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor',[0.5,0.5,0.5]);
axis square; box off;
h=get(gcf, 'currentaxes');
set(h, 'fontsize', 15, 'linewidth', 2);
box off;
get(gca,'ticklength');
set(gca,'TickLength',[0.03,0.02]);
set(gca,'TickDir','out');
ylim([0 max([rate_in_without0,rate_out_without0])]);
xlabel('path'); ylabel('frequency(Hz)');

subplot(3,3,6);
plot(rate_out_without0,'k-o','LineWidth',2,'MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor',[0.5,0.5,0.5]);
axis square; box off;
h=get(gcf, 'currentaxes');
set(h, 'fontsize', 15, 'linewidth', 2);
box off;
get(gca,'ticklength');
set(gca,'TickLength',[0.03,0.02]);
set(gca,'TickDir','out');
ylim([0 max([rate_in_without0,rate_out_without0])]);
xlabel('path'); ylabel('frequency(Hz)');

subplot(3,3,[8 9]);
scatter(idx_left,rate(compartment_side == 10),150,'or','filled');
hold on
scatter(idx_right,rate(compartment_side == -10),150,'ob','filled');
scatter(idx_door,rate(compartment_side == 1),150,'og','filled');
box off;
h=get(gcf, 'currentaxes');
set(h, 'fontsize', 15, 'linewidth', 2);
box off;
get(gca,'ticklength');
set(gca,'TickLength',[0.03,0.02]);
set(gca,'TickDir','out');
ylim([0 max(rate)]);
xlabel('path'); ylabel('frequency(Hz)');
legend({'left','right','door'},'Position',[.1 .1 .2 .2]);
%%

set(temporal_fig,'units','normalized','outerposition',[0 0 1 1 ]);
set(temporal_fig,'Units','Inches');
pos = get(temporal_fig,'Position');
set(temporal_fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])


    if save_figure == 1
            complete_name = strcat(pathfigure,inputfile(1:end-4),'_t',num2str(p),'_c',num2str(i),'_temporal_activity_two_comp');
            print('-dpng',complete_name,'-r100');

   end
   if close_figure == 1
        close(temporal_fig);
   end


%% Save results
results_txt_2comp(results_folder,inputfile,p,i,peak_rate_global,peak_rate_left,peak_rate_right,...
    peak_rate_door,corrcoeff,rotated_corrcoeff,PFD_deg,Vector_length,peak_HD,PFD_deg_left,...
    Vector_length_left,peak_HD_left,PFD_deg_right,Vector_length_right,peak_HD_right,...
    PFD_deg_door,Vector_length_door,peak_HD_door,Vector_length_double,flipping_score,...
    flipping_score_range,Vector_length_quadri,Four_dir_score,Four_dir_score_minmax,...
    Four_dir_score_range,Four_dir_score_minmax_range,ret,peak_frequency_theta,theta_power,...
    delta_power,sig_threshold,theta_modulation,number_of_corr_peaks);
%% save mat
nameformat = strcat(mat_folder,filename(1:end-4),'_t',num2str(p),'_c',num2str(i),'.mat');
save(nameformat,...
    'TTMtx2','PosMtx','filename','ts','xy_left','xy_right','xy_door','left','right','door',...
    'waveforms','spike_left','spike_right','spike_door','ratemap_left','ratemap_right',...
    'ratemap_comp1','rotated_ratemap_comp1','ratemap_comp2',...
    'dirrate','dirrate_left','dirrate_right','dirrate_door','corrcoeff','rotated_corrcoeff');

%% Splitting 4 compartments
    elseif compartment == 4

% [spkx_in1,spky_in1,spkx_in2,spky_in2,spkx_in3,spky_in3,spkx_in4,spky_in4,spkx,spky,ts_in1,ts_in2,ts_in3,ts_in4,posgx_in1,posgy_in1,posgx_in2,posgy_in2,posgx_in3,posgy_in3,posgx_in4,posgy_in4,rate_in1,rate_in2,rate_in3,rate_in4,ratemap_in1,ratemap_in2,ratemap_in3,ratemap_in4,rotated_ratemap_in1,rotated_ratemap_in2,rotated_ratemap_in3,rotated_ratemap_in4,RrankO1,RrankO2,RrankO3,RrankO4,RrankO5,RrankO6,rotated_RrankO1,rotated_RrankO2,rotated_RrankO3,rotated_RrankO4,rotated_RrankO5,rotated_RrankO6] = Splitting_four_compartment_activities(ts,posts,posrx,posry,posgx,posgy,posrx_in1,posry_in1,posrx_in2,posry_in2,posrx_in3,posry_in3,posrx_in4,posry_in4,x,x2,x3,x4,y,y2,y3,y4,pos_in_rectangle1,pos_in_rectangle2,pos_in_rectangle3,pos_in_rectangle4,posts_in1,posts_in2,posts_in3,posts_in4,binWidth,shape,V,V2,V3,V4,xx,xx2,xx3,xx4);
[spkx,spky,spike_top_left,spike_bottom_left,spike_top_right,spike_bottom_right,door_top,door_left,door_bottom,door_right] = Splitting_four_compartments_and_doors(ts,posx,posy,posts,xy_top_left,xy_bottom_left,xy_top_right,xy_bottom_right,xy_top_door,xy_left_door,xy_bottom_door,xy_right_door);

%% open figure
four_comp = figure('visible', 'off');
%% Rate map global
[posx,posy] = smooth_path(posx,posy); [posx,posy] = center_path(posx,posy,shape);
[mapAxis] = mapaxis(posx,posy,binWidth); visited = visitedBins(posx,posy,mapAxis);   
[spkx,spky] = get_pos_spikes(ts,posx,posy,posts); % get the position of the spikes

subplot(4,4,1);
plot(posx,posy,'color',[.2 .2 .2],'linewidth',2);
hold on
scatter(spkx,spky,200,'.r');
set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);
axis off

[ratemap] = ratemap_gaussian(2*binWidth,spkx,spky,posx,posy,posts,binWidth,mapAxis); % build matrix ratemap 
ratemap(visited==0) = NaN;
subplot(4,4,2);
myfig=pcolor(ratemap);
colormap jet
set(myfig,'Edgecolor','none')
axis xy
set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);
axis off
box off
peak_rate=max(ratemap(:))*1000; peak_rate_global=peak_rate;peak_rate_global=sprintf('%.2f', peak_rate_global);
title(strcat('Peak rate = ',peak_rate_global, ' Hz'));
peak_rate_global = peak_rate;

%% polar plot global
posts = PosMtx(:,1);
posgx = PosMtx(:,2);
posgy = PosMtx(:,3);
[posgx,posgy] = smooth_path(posgx,posgy);
posts = PosMtx(:,1);
posrx = PosMtx(:,4);
posry = PosMtx(:,5);
[posrx,posry] = smooth_path(posrx,posry);

subplot(4,4,3);
[PFD_deg,Vector_length,peak_HD,dirrate,deg_theta,posdirts,dirmap]= polarplot_compartment_new(posts,posgx,posgy,posrx,posry,ts,correction,binSizeDir);
peak_HD_2=sprintf('%.2f', peak_HD);
Vector_length_2=sprintf('%.2f', Vector_length);
set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);
title({
    ['Pic = ', peak_HD_2,'Hz, Vector length = ',Vector_length_2]},'fontsize',12);
set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);xlabel('global')

%path for each compartment and the door
subplot(4,4,4);
plot(posrx,posry);
hold on
% plot position of position in each quadrant with different color
scatter(posrx(topleft),posry(topleft),200,'r.'); % points inside
scatter(posrx(bottomleft),posry(bottomleft),200,'g.'); % points inside
scatter(posrx(bottomright),posry(bottomright),200,'b.'); 
scatter(posrx(topright),posry(topright),200,'y.'); 
scatter(posrx(topdoor),posry(topdoor),200,'.','MarkerEdgeColor',[1 0 1]); % points inside
scatter(posrx(leftdoor),posry(leftdoor),200,'.','MarkerEdgeColor',[.3 .3 .3]); % points inside
scatter(posrx(bottomdoor),posry(bottomdoor),200,'.','MarkerEdgeColor',[1 .5 0]); 
scatter(posrx(rightdoor),posry(rightdoor),200,'.','MarkerEdgeColor',[.5 0 0]); 
axis square
axis off

subplot(4,4,13);
% polar plot left compartment
[PFD_deg_top_left,Vector_length_top_left,peak_HD_top_left,dirrate_top_left,deg_theta,posdirts,dirmap]= polarplot_compartment_new(posts(topleft),posgx(topleft),posgy(topleft),posrx(topleft),posry(topleft),ts(spike_top_left),correction,binSizeDir);
peak_HD_2=sprintf('%.2f', peak_HD_top_left);
Vector_length_2=sprintf('%.2f', Vector_length_top_left);
set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);
title({['Pic = ',peak_HD_2,'Hz, Vector length=',Vector_length_2]});
set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);

subplot(4,4,14);
% polar plot right compartment
[PFD_deg_bottom_left,Vector_length_bottom_left,peak_HD_bottom_left,dirrate_bottom_left,deg_theta,posdirts,dirmap]= polarplot_compartment_new(posts(bottomleft),posgx(bottomleft),posgy(bottomleft),posrx(bottomleft),posry(bottomleft),ts(spike_bottom_left),correction,binSizeDir);
peak_HD_2=sprintf('%.2f', peak_HD_bottom_left);
Vector_length_2=sprintf('%.2f', Vector_length_bottom_left);
set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);
title({['Pic = ',peak_HD_2,'Hz, Vector length=',Vector_length_2]});
set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);

subplot(4,4,15);
% polar plot door
[PFD_deg_bottom_right,Vector_length_bottom_right,peak_HD_bottom_right,dirrate_bottom_right,deg_theta,posdirts,dirmap]= polarplot_compartment_new(posts(bottomright),posgx(bottomright),posgy(bottomright),posrx(bottomright),posry(bottomright),ts(spike_bottom_right),correction,binSizeDir);
peak_HD_2=sprintf('%.2f', peak_HD_bottom_right);
Vector_length_2=sprintf('%.2f', Vector_length_bottom_right);
set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);
title({['Pic = ',peak_HD_2,'Hz, Vector length=',Vector_length_2]});
set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);

subplot(4,4,16);
% polar plot door
[PFD_deg_top_right,Vector_length_top_right,peak_HD_top_right,dirrate_top_right,deg_theta,posdirts,dirmap]= polarplot_compartment_new(posts(topright),posgx(topright),posgy(topright),posrx(topright),posry(topright),ts(spike_top_right),correction,binSizeDir);
peak_HD_2=sprintf('%.2f', peak_HD_top_right);
Vector_length_2=sprintf('%.2f', Vector_length_top_right);
set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);
title({['Pic = ',peak_HD_2,'Hz, Vector length=',Vector_length_2]});
set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);

% ratemap compartments
[posrx,posry] = center_path(posrx,posry,shape);

subplot(4,4,5);
plot(posrx,posry,'color',[.2 .2 .2],'linewidth',2);
[mapAxis] = mapaxis(posrx,posry,binWidth); visited = visitedBins(posrx(topleft),posry(topleft),mapAxis);
ts_top_left = ts(spike_top_left);
[spkx_top_left,spky_top_left] = get_pos_spikes(ts_top_left,posrx(topleft),posry(topleft),posts(topleft)); % get the position of the spikes
hold on
scatter(spkx_top_left,spky_top_left,200,'.r');
set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);
axis off

[ratemap_top_left] = ratemap_gaussian(2*binWidth,spkx_top_left,spky_top_left,posrx(topleft),posry(topleft),posts(topleft),binWidth,mapAxis); % build matrix ratemap 
ratemap_top_left(visited==0) = NaN;
subplot(4,4,9);
myfig=pcolor(ratemap_top_left);
colormap jet
set(myfig,'Edgecolor','none')
axis xy
set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);
axis off
box off
peak_rate=max(ratemap_top_left(:))*1000; peak_ratetop_left=peak_rate;peak_ratetop_left=sprintf('%.2f', peak_ratetop_left);
title(strcat('Peak rate = ',peak_ratetop_left, ' Hz'));
peak_rate_top_left = peak_rate;

subplot(4,4,6);
plot(posrx,posry,'color',[.2 .2 .2],'linewidth',2);
[mapAxis] = mapaxis(posrx,posry,binWidth); visited = visitedBins(posrx(bottomleft),posry(bottomleft),mapAxis);
ts_bottom_left = ts(spike_bottom_left);
[spkx_bottom_left,spky_bottom_left] = get_pos_spikes(ts_bottom_left,posrx(bottomleft),posry(bottomleft),posts(bottomleft)); % get the position of the spikes
hold on
scatter(spkx_bottom_left,spky_bottom_left,200,'.r');
set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);
axis off

[ratemap_bottom_left] = ratemap_gaussian(2*binWidth,spkx_bottom_left,spky_bottom_left,posrx(bottomleft),posry(bottomleft),posts(bottomleft),binWidth,mapAxis); % build matrix ratemap 
ratemap_bottom_left(visited==0) = NaN;
subplot(4,4,10);
myfig=pcolor(ratemap_bottom_left);
colormap jet
set(myfig,'Edgecolor','none')
axis xy
set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);
axis off
box off
peak_rate=max(ratemap_bottom_left(:))*1000; peak_rate_bottom_left=peak_rate;peak_rate_bottom_left=sprintf('%.2f', peak_rate_bottom_left);
title(strcat('Peak rate = ',peak_rate_bottom_left, ' Hz'));
peak_rate_bottom_left = peak_rate;


subplot(4,4,7);
plot(posrx,posry,'color',[.2 .2 .2],'linewidth',2);
[mapAxis] = mapaxis(posrx,posry,binWidth); visited = visitedBins(posrx(bottomright),posry(bottomright),mapAxis);
ts_bottom_right = ts(spike_bottom_right);
[spkx_bottom_right,spky_bottom_right] = get_pos_spikes(ts_bottom_right,posrx(bottomright),posry(bottomright),posts(bottomright)); % get the position of the spikes
hold on
scatter(spkx_bottom_right,spky_bottom_right,200,'.r');
set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);
axis off

[ratemap_bottom_right] = ratemap_gaussian(2*binWidth,spkx_bottom_right,spky_bottom_right,posrx(bottomright),posry(bottomright),posts(bottomright),binWidth,mapAxis); % build matrix ratemap 
ratemap_bottom_right(visited==0) = NaN;
subplot(4,4,11);
myfig=pcolor(ratemap_bottom_right);
colormap jet
set(myfig,'Edgecolor','none')
axis xy
set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);
axis off
box off
peak_rate=max(ratemap_bottom_right(:))*1000; peak_rate_bottom_right=peak_rate;peak_rate_bottom_right=sprintf('%.2f', peak_rate_bottom_right);
title(strcat('Peak rate = ',peak_rate_bottom_right, ' Hz'));
peak_rate_bottom_right = peak_rate;


subplot(4,4,8);
plot(posrx,posry,'color',[.2 .2 .2],'linewidth',2);
[mapAxis] = mapaxis(posrx,posry,binWidth); visited = visitedBins(posrx(topright),posry(topright),mapAxis);
ts_top_right = ts(spike_top_right);
[spkx_top_right,spky_top_right] = get_pos_spikes(ts_top_right,posrx(topright),posry(topright),posts(topright)); % get the position of the spikes
hold on
scatter(spkx_top_right,spky_top_right,200,'.r');
set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);
axis off

[ratemap_top_right] = ratemap_gaussian(2*binWidth,spkx_top_right,spky_top_right,posrx(topright),posry(topright),posts(topright),binWidth,mapAxis); % build matrix ratemap 
ratemap_top_right(visited==0) = NaN;
subplot(4,4,12);
myfig=pcolor(ratemap_top_right);
colormap jet
set(myfig,'Edgecolor','none')
axis xy
set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);
axis off
box off
peak_rate=max(ratemap_top_right(:))*1000; peak_rate_top_right=peak_rate;peak_rate_top_right=sprintf('%.2f', peak_rate_top_right);
title(strcat('Peak rate = ',peak_rate_top_right, ' Hz'));
peak_rate_top_right = peak_rate;


%% save figures of 4 compartments polar plots  
set(four_comp,'units','normalized','outerposition',[0 0 1 1 ]);
set(four_comp,'Units','Inches');
pos = get(four_comp,'Position');
set(four_comp,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])


    if save_figure == 1
            complete_name = strcat(pathfigure,inputfile(1:end-4),'_t',num2str(p),'_c',num2str(i),'_split_four_comp');
            print('-dpng',complete_name,'-r100');

   end
   if close_figure == 1
        close(four_comp);
   end

%%
general_fig_4comp = figure('visible', 'off');
[selected_el,average_selected_el,std_plus_selected_el,std__selected_el,std_minus_selected_el,...
    average_derive_selected_el,std__derive_selected_el,std_plus_derive_selected_el,std_minus_derive_selected_el]= waveform(waveforms);
subplot(2,3,1);
    plot(average_selected_el,'linewidth',3,'color','k'); axis square;
    hold on
    plot(std_plus_selected_el,'--','linewidth',3,'color','k'); axis square;
    plot(std_minus_selected_el,'--','linewidth',3,'color','k'); axis square;
    hold off
    axis off
    box off
    title('Mean waveform','fontsize',12);
    xlim([0 32]);
    set(gca, 'XTickLabel', [])
    set(gca, 'XTick', [])
    set(gca, 'YTickLabel', [])
    set(gca, 'YTick', [])
    h=get(gcf, 'currentaxes');
    set(h, 'fontsize', 12, 'linewidth', 2);

subplot(2,3,2);
plot(posx,posy,'color',[.4 .4 .4],'linewidth',2);%,spkx,spky,'.r');
 hold on
scatter(spkx,spky,400,'.r');
set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);
axis off
box off

posts = PosMtx(:,1);
posgx = PosMtx(:,2);
posgy = PosMtx(:,3);
[posgx,posgy] = smooth_path(posgx,posgy);
posts = PosMtx(:,1);
posrx = PosMtx(:,4);
posry = PosMtx(:,5);
[posrx,posry] = smooth_path(posrx,posry);
subplot(2,3,3);
polarplot_compartment_new(posts,posgx,posgy,posrx,posry,ts,correction,binSizeDir);
set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);

% calculate and plot cell firing/speed relationship
subplot(2,3,4);
[ret,beta,f0]=speed_firing(posx,posy,posts,ts);
h=get(gcf, 'currentaxes');
set(h, 'fontsize', 15, 'linewidth', 2);
box off;
get(gca,'ticklength');
set(gca,'TickLength',[0.03,0.02]);
set(gca,'TickDir','out');
    
% calculate and plot time firing autocorrelagram
subplot(2,3,5);
[autocorr,y] = autocorrelogram(TimeWindow,Tbin,ts);
h=get(gcf, 'currentaxes');
set(h, 'fontsize', 15, 'linewidth', 2);
box off;
get(gca,'ticklength');
set(gca,'TickLength',[0.03,0.02]);
set(gca,'TickDir','out');
title('Autocorrelogram','fontsize',12);

subplot(2,3,6);
[mean_freq,freq_peak_distr,peak_theta,peak_delta] = intrinsicfiring(TimeWindow,Tbin,autocorr,y);
title(strcat('Peak theta=',num2str(sprintf('%.2f', freq_peak_distr),'Hz')),'fontsize',12);
[peak_frequency_theta,theta_power,delta_power,sig_threshold,theta_modulation] = theta_delta(TimeWindow,Tbin,ts);
set(gca,'TickDir','out','Box','off');
h=get(gcf, 'currentaxes');
set(h, 'fontsize', 12, 'linewidth', 2);
set(gca,'TickLength',[0.02,0.01]);
title(strcat('Peak theta=',num2str(sprintf('%.2f', freq_peak_distr),'Hz')),'fontsize',12);

%% save general figures 
set(general_fig_4comp,'units','normalized','outerposition',[0 0 1 1 ]);
set(general_fig_4comp,'Units','Inches');
pos = get(general_fig_4comp,'Position');
set(general_fig_4comp,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

   if save_figure == 1
            complete_name = strcat(pathfigure,inputfile(1:end-4),'_t',num2str(p),'_c',num2str(i),'_general_four_comp');
            print('-dpng',complete_name,'-r100');
%             print(two_comp,complete_name,'-dpdf')
   end
    if close_figure == 1
        close(general_fig_4comp);
    end



%% open figure
four_door = figure('visible', 'off');
%% Rate map global
[posx,posy] = smooth_path(posx,posy); [posx,posy] = center_path(posx,posy,shape);
[mapAxis] = mapaxis(posx,posy,binWidth); visited = visitedBins(posx,posy,mapAxis);   
[spkx,spky] = get_pos_spikes(ts,posx,posy,posts); % get the position of the spikes

subplot(4,4,1);
plot(posx,posy,'color',[.2 .2 .2],'linewidth',2);
hold on
scatter(spkx,spky,200,'.r');
set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);
axis off

[ratemap] = ratemap_gaussian(2*binWidth,spkx,spky,posx,posy,posts,binWidth,mapAxis); % build matrix ratemap 
ratemap(visited==0) = NaN;
subplot(4,4,2);
myfig=pcolor(ratemap);
colormap jet
set(myfig,'Edgecolor','none')
axis xy
set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);
axis off
box off
peak_rate=max(ratemap(:))*1000; peak_rate2=peak_rate;peak_rate2=sprintf('%.2f', peak_rate2);
title(strcat('Peak rate = ',peak_rate2, ' Hz'));

%% polar plot global
posts = PosMtx(:,1);
posgx = PosMtx(:,2);
posgy = PosMtx(:,3);
[posgx,posgy] = smooth_path(posgx,posgy);
posts = PosMtx(:,1);
posrx = PosMtx(:,4);
posry = PosMtx(:,5);
[posrx,posry] = smooth_path(posrx,posry);

subplot(4,4,3);
polarplot_compartment_new(posts,posgx,posgy,posrx,posry,ts,correction,binSizeDir);
peak_HD_2=sprintf('%.2f', peak_HD);
Vector_length_2=sprintf('%.2f', Vector_length);
set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);
title({
    ['Pic = ', peak_HD_2,'Hz, Vector length = ',Vector_length_2]},'fontsize',12);
set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);xlabel('global')

%path for each compartment and the door
subplot(4,4,4);
plot(posrx,posry);
hold on
% plot position of position in each quadrant with different color
scatter(posrx(topleft),posry(topleft),200,'r.'); % points inside
scatter(posrx(bottomleft),posry(bottomleft),200,'g.'); % points inside
scatter(posrx(bottomright),posry(bottomright),200,'b.'); 
scatter(posrx(topright),posry(topright),200,'y.'); 
scatter(posrx(topdoor),posry(topdoor),200,'.','MarkerEdgeColor',[1 0 1]); % points inside
scatter(posrx(leftdoor),posry(leftdoor),200,'.','MarkerEdgeColor',[.3 .3 .3]); % points inside
scatter(posrx(bottomdoor),posry(bottomdoor),200,'.','MarkerEdgeColor',[1 .5 0]); 
scatter(posrx(rightdoor),posry(rightdoor),200,'.','MarkerEdgeColor',[.5 0 0]); 
axis square
axis off

subplot(4,4,13);
% polar plot top door
[PFD_deg_door_top,Vector_length_door_top,peak_HD_door_top,dirrate_door_top,deg_theta,posdirts,dirmap]= polarplot_compartment_new(posts(topdoor),posgx(topdoor),posgy(topdoor),posrx(topdoor),posry(topdoor),ts(door_top),correction,binSizeDir);
peak_HD_2=sprintf('%.2f', peak_HD_door_top);
Vector_length_2=sprintf('%.2f', Vector_length_door_top);
set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);
title({['Pic = ',peak_HD_2,'Hz, Vector length=',Vector_length_2]});
set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);

subplot(4,4,14);
% polar plot left door
[PFD_deg_door_left,Vector_length_door_left,peak_HD_door_left,dirrate_door_left,deg_theta,posdirts,dirmap]= polarplot_compartment_new(posts(leftdoor),posgx(leftdoor),posgy(leftdoor),posrx(leftdoor),posry(leftdoor),ts(door_left),correction,binSizeDir);
peak_HD_2=sprintf('%.2f', peak_HD_door_left);
Vector_length_2=sprintf('%.2f', Vector_length_door_left);
set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);
title({['Pic = ',peak_HD_2,'Hz, Vector length=',Vector_length_2]});
set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);

subplot(4,4,15);
% polar plot bottom door
[PFD_deg_door_bottom,Vector_length_door_bottom,peak_HD_door_bottom,dirrate_door_bottom,deg_theta,posdirts,dirmap]= polarplot_compartment_new(posts(bottomdoor),posgx(bottomdoor),posgy(bottomdoor),posrx(bottomdoor),posry(bottomdoor),ts(door_bottom),correction,binSizeDir);
peak_HD_2=sprintf('%.2f', peak_HD_door_bottom);
Vector_length_2=sprintf('%.2f', Vector_length_door_bottom);
set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);
title({['Pic = ',peak_HD_2,'Hz, Vector length=',Vector_length_2]});
set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);

subplot(4,4,16);
% polar plot right door
[PFD_deg_door_right,Vector_length_door_right,peak_HD_door_right,dirrate_door_right,deg_theta,posdirts,dirmap]= polarplot_compartment_new(posts(rightdoor),posgx(rightdoor),posgy(rightdoor),posrx(rightdoor),posry(rightdoor),ts(door_right),correction,binSizeDir);
peak_HD_2=sprintf('%.2f', peak_HD_door_right);
Vector_length_2=sprintf('%.2f', Vector_length_door_right);
set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);
title({['Pic = ',peak_HD_2,'Hz, Vector length=',Vector_length_2]});
set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);

% ratemap door
[posrx,posry] = center_path(posrx,posry,shape);

subplot(4,4,5);
plot(posrx,posry,'color',[.2 .2 .2],'linewidth',2);
[mapAxis] = mapaxis(posrx,posry,binWidth); visited = visitedBins(posrx(topdoor),posry(topdoor),mapAxis);
ts_top_left = ts(door_left);
[spkx_top_left,spky_top_left] = get_pos_spikes(ts_top_left,posrx(topdoor),posry(topdoor),posts(topdoor)); % get the position of the spikes
hold on
scatter(spkx_top_left,spky_top_left,200,'.r');
set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);
axis off

[ratemap_door_top] = ratemap_gaussian(2*binWidth,spkx_top_left,spky_top_left,posrx(topdoor),posry(topdoor),posts(topdoor),binWidth,mapAxis); % build matrix ratemap 
ratemap_door_top(visited==0) = NaN;
subplot(4,4,9);
myfig=pcolor(ratemap_door_top);
colormap jet
set(myfig,'Edgecolor','none')
axis xy
set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);
axis off
box off
peak_rate=max(ratemap_door_top(:))*1000; peak_rate2=peak_rate;peak_rate2=sprintf('%.2f', peak_rate2);
title(strcat('Peak rate = ',peak_rate2, ' Hz'));
peak_rate_door_top = peak_rate;

subplot(4,4,6);
plot(posrx,posry,'color',[.2 .2 .2],'linewidth',2);
[mapAxis] = mapaxis(posrx,posry,binWidth); visited = visitedBins(posrx(leftdoor),posry(leftdoor),mapAxis);
ts_bottom_left = ts(door_left);
[spkx_bottom_left,spky_bottom_left] = get_pos_spikes(ts_bottom_left,posrx(leftdoor),posry(leftdoor),posts(leftdoor)); % get the position of the spikes
hold on
scatter(spkx_bottom_left,spky_bottom_left,200,'.r');
set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);
axis off

[ratemap_door_left] = ratemap_gaussian(2*binWidth,spkx_bottom_left,spky_bottom_left,posrx(leftdoor),posry(leftdoor),posts(leftdoor),binWidth,mapAxis); % build matrix ratemap 
ratemap_door_left(visited==0) = NaN;
subplot(4,4,10);
myfig=pcolor(ratemap_door_left);
colormap jet
set(myfig,'Edgecolor','none')
axis xy
set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);
axis off
box off
peak_rate=max(ratemap_door_left(:))*1000; peak_rate2=peak_rate;peak_rate2=sprintf('%.2f', peak_rate2);
title(strcat('Peak rate = ',peak_rate2, ' Hz'));
peak_rate_door_left = peak_rate;

subplot(4,4,7);
plot(posrx,posry,'color',[.2 .2 .2],'linewidth',2);
[mapAxis] = mapaxis(posrx,posry,binWidth); visited = visitedBins(posrx(bottomdoor),posry(bottomdoor),mapAxis);
ts_bottom_right = ts(door_bottom);
[spkx_bottom_right,spky_bottom_right] = get_pos_spikes(ts_bottom_right,posrx(bottomdoor),posry(bottomdoor),posts(bottomdoor)); % get the position of the spikes
hold on
scatter(spkx_bottom_right,spky_bottom_right,200,'.r');
set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);
axis off

[ratemap_door_bottom] = ratemap_gaussian(2*binWidth,spkx_bottom_right,spky_bottom_right,posrx(bottomdoor),posry(bottomdoor),posts(bottomdoor),binWidth,mapAxis); % build matrix ratemap 
ratemap_door_bottom(visited==0) = NaN;
subplot(4,4,11);
myfig=pcolor(ratemap_door_bottom);
colormap jet;
set(myfig,'Edgecolor','none')
axis xy
set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);
axis off
box off
peak_rate=max(ratemap_door_bottom(:))*1000; peak_rate2=peak_rate;peak_rate2=sprintf('%.2f', peak_rate2);
title(strcat('Peak rate = ',peak_rate2, ' Hz'));
peak_rate_door_bottom = peak_rate;

subplot(4,4,8);
plot(posrx,posry,'color',[.2 .2 .2],'linewidth',2);
[mapAxis] = mapaxis(posrx,posry,binWidth); visited = visitedBins(posrx(rightdoor),posry(rightdoor),mapAxis);
ts_top_right = ts(door_right);
[spkx_top_right,spky_top_right] = get_pos_spikes(ts_top_right,posrx(rightdoor),posry(rightdoor),posts(rightdoor)); % get the position of the spikes
hold on
scatter(spkx_top_right,spky_top_right,200,'.r');
set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);
axis off

[ratemap_door_right] = ratemap_gaussian(2*binWidth,spkx_top_right,spky_top_right,posrx(rightdoor),posry(rightdoor),posts(rightdoor),binWidth,mapAxis); % build matrix ratemap 
ratemap_door_right(visited==0) = NaN;
subplot(4,4,12);
myfig=pcolor(ratemap_door_right);
colormap jet;
set(myfig,'Edgecolor','none')
axis xy;
set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);
axis off;
box off;
peak_rate=max(ratemap_door_right(:))*1000; peak_rate2=peak_rate;peak_rate2=sprintf('%.2f', peak_rate2);
title(strcat('Peak rate = ',peak_rate2, ' Hz'));
peak_rate_door_right = peak_rate;


%% save figures of 4 compartments polar plots  
set(four_door,'units','normalized','outerposition',[0 0 1 1 ]);
set(four_door,'Units','Inches');
pos = get(four_door,'Position');
set(four_door,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

 if save_figure == 1
            complete_name = strcat(pathfigure,inputfile(1:end-4),'_t',num2str(p),'_c',num2str(i),'_four_comp');
            print('-dpng',complete_name,'-r100');

   end
   if close_figure == 1
        close(four_door);
   end
%%

colored_plot = figure('visible', 'off');
[dirrate] = polarplot_v2_colored_new(PosMtx,spkx,spky,ts,binSizeDir,correction,shape);

[Vector_length_double,flipping_score,flipping_score_range,Vector_length_quadri,Four_dir_score,Four_dir_score_minmax,Four_dir_score_range,Four_dir_score_minmax_range,correlation,minimums,maximums,corr_value,corr_bin]=flipping_scores(dirrate,binSizeDir,PFD_deg);
minimums_angles = corr_bin(minimums);
maximums_angles = corr_bin(maximums);
minimums_values = corr_value(minimums);
maximums_values = corr_value(maximums);

number_of_corr_peaks = sum(minimums)+sum(maximums);

subplot(2,2,3);
plot(correlation(:,2),correlation(:,1),'k','linewidth',2); hold on;
scatter(corr_bin(maximums),corr_value(maximums),100,'o','markeredgecolor',[.3 .3 .3],'markerfacecolor',[.7 .7 .7]);
scatter(corr_bin(minimums),corr_value(minimums),100,'o','markeredgecolor',[.3 .3 .3],'markerfacecolor',[.7 .7 .7]);
axis square; box off;
h=get(gcf, 'currentaxes');
set(h, 'fontsize', 15, 'linewidth', 2);
box off;
get(gca,'ticklength');
set(gca,'TickLength',[0.03,0.02]);
set(gca,'TickDir','out');
xlim([0 360]);
xticks(0:45:360);
xlabel('degre'); ylabel('correlation coefficient');

subplot(2,2,4); axis off;
text(0,1,strcat('flip score = ',num2str(flipping_score)), 'fontsize', 20);
text(0,0.8,strcat('flips score range = ',num2str(flipping_score_range)), 'fontsize', 20);
text(0,0.6,strcat('quadiscore = ',num2str(Four_dir_score)), 'fontsize', 20);
text(0,0.4,strcat('quadiscore range = ',num2str(Four_dir_score_range)), 'fontsize', 20);
text(0,0.2,strcat('quadiscore minmax = ',num2str(Four_dir_score_minmax)), 'fontsize', 20);
text(0,0,strcat('quadiscore minmax range = ',num2str(Four_dir_score_minmax_range)), 'fontsize', 20);


set(colored_plot,'units','normalized','outerposition',[0 0 1 1 ]);
set(colored_plot,'Units','Inches');
pos = get(colored_plot,'Position');
set(colored_plot,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

   if save_figure == 1
            complete_name = strcat(pathfigure,inputfile(1:end-4),'_t',num2str(p),'_c',num2str(i),'_colored_plot_four_comp');
            print('-dpng',complete_name,'-r100');
%             print(colored_plot,complete_name,'-dpdf');

   end
   if close_figure == 1
        close(colored_plot);
   end
%% spatial correlation between compartments non rotated and rotated
% ratemp comp 1
[posrx_comp1,posry_comp1] = center_path(posrx(topleft),posry(topleft),shape);
[mapAxis_comp1] = mapaxis(posrx_comp1,posry_comp1,binWidth);
visited_comp1 = visitedBins(posrx_comp1,posry_comp1,mapAxis_comp1); % calculate the visited zones of the arena 
[spkx_comp1,spky_comp1] = get_pos_spikes(ts(spike_top_left),posrx_comp1,posry_comp1,posts(topleft)); % get the position of the spikes
[ratemap_comp1] = ratemap_gaussian(2*binWidth,spkx_comp1,spky_comp1,posrx_comp1,posry_comp1,posts(topleft),binWidth,mapAxis_comp1); % build matrix ratemap 
ratemap_comp1(visited_comp1==0) = NaN;
ratemap_top_left2 = ratemap_comp1;

%ratemp comp 2
[posrx_comp2,posry_comp2] = center_path(posrx(bottomleft),posry(bottomleft),shape);
[mapAxis_comp2] = mapaxis(posrx_comp2,posry_comp2,binWidth);
visited_comp2 = visitedBins(posrx_comp2,posry_comp2,mapAxis_comp2); % calculate the visited zones of the arena 
[spkx_comp2,spky_comp2] = get_pos_spikes(ts(spike_bottom_left),posrx_comp2,posry_comp2,posts(bottomleft)); % get the position of the spikes
[ratemap_comp2] = ratemap_gaussian(2*binWidth,spkx_comp2,spky_comp2,posrx_comp2,posry_comp2,posts(bottomleft),binWidth,mapAxis_comp2); % build matrix ratemap 
ratemap_comp2(visited_comp2==0) = NaN;
ratemap_down_left2 = ratemap_comp2;

%ratemp comp 3
[posrx_comp3,posry_comp3] = center_path(posrx(bottomright),posry(bottomright),shape);
[mapAxis_comp3] = mapaxis(posrx_comp3,posry_comp3,binWidth);
visited_comp3 = visitedBins(posrx_comp3,posry_comp3,mapAxis_comp3); % calculate the visited zones of the arena 
[spkx_comp3,spky_comp3] = get_pos_spikes(ts(spike_bottom_right),posrx_comp3,posry_comp3,posts(bottomright)); % get the position of the spikes
[ratemap_comp3] = ratemap_gaussian(2*binWidth,spkx_comp3,spky_comp3,posrx_comp3,posry_comp3,posts(bottomright),binWidth,mapAxis_comp3); % build matrix ratemap 
ratemap_comp3(visited_comp3==0) = NaN;
ratemap_down_right2 = ratemap_comp3;

%ratemp comp 4
[posrx_comp4,posry_comp4] = center_path(posrx(topright),posry(topright),shape);
[mapAxis_comp4] = mapaxis(posrx_comp4,posry_comp4,binWidth);
visited_comp4 = visitedBins(posrx_comp4,posry_comp4,mapAxis_comp4); % calculate the visited zones of the arena 
[spkx_comp4,spky_comp4] = get_pos_spikes(ts(spike_top_right),posrx_comp4,posry_comp4,posts(topright)); % get the position of the spikes
[ratemap_comp4] = ratemap_gaussian(2*binWidth,spkx_comp4,spky_comp4,posrx_comp4,posry_comp4,posts(topright),binWidth,mapAxis_comp4); % build matrix ratemap 
ratemap_comp4(visited_comp4==0) = NaN;
ratemap_top_right2 = ratemap_comp4;

minLength = min([length(ratemap_top_left2), length(ratemap_down_left2), length(ratemap_down_right2), length(ratemap_top_right2)]);
ratemap_top_left2 = ratemap_top_left2(1:minLength,1:minLength); ratemap_down_left2=ratemap_down_left2(1:minLength,1:minLength); ratemap_down_right2=ratemap_down_right2(1:minLength,1:minLength); ratemap_top_right2=ratemap_top_right2(1:minLength,1:minLength); 

[corrcoeff1,nonZero_corrcoeff1,RrankO1] = get_2D_correlation_value(ratemap_top_left2,ratemap_down_left2);
[corrcoeff2,nonZero_corrcoeff2,RrankO2] = get_2D_correlation_value(ratemap_top_left2,ratemap_down_right2);
[corrcoeff3,nonZero_corrcoeff3,RrankO3] = get_2D_correlation_value(ratemap_top_left2,ratemap_top_right2);
[corrcoeff4,nonZero_corrcoeff4,RrankO4] = get_2D_correlation_value(ratemap_down_left2,ratemap_down_right2);
[corrcoeff5,nonZero_corrcoeff5,RrankO5] = get_2D_correlation_value(ratemap_down_left2,ratemap_top_right2);
[corrcoeff6,nonZero_corrcoeff6,RrankO6] = get_2D_correlation_value(ratemap_down_right2,ratemap_top_right2);

meancorrcoef = mean([corrcoeff1,corrcoeff2,corrcoeff3,corrcoeff4,corrcoeff5,corrcoeff6]);

%ratemp comp 2 rotated 90°
[posrx_comp2,posry_comp2] = center_path(posry(bottomleft),posrx(bottomleft),shape);
posry_comp2=-posry_comp2;
[mapAxis_comp2] = mapaxis(posrx_comp2,posry_comp2,binWidth);
visited_comp2 = visitedBins(posrx_comp2,posry_comp2,mapAxis_comp2); % calculate the visited zones of the arena 
[spkx_comp2,spky_comp2] = get_pos_spikes(ts(spike_bottom_left),posrx_comp2,posry_comp2,posts(bottomleft)); % get the position of the spikes
[ratemap_comp2] = ratemap_gaussian(2*binWidth,spkx_comp2,spky_comp2,posrx_comp2,posry_comp2,posts(bottomleft),binWidth,mapAxis_comp2); % build matrix ratemap 
ratemap_comp2(visited_comp2==0) = NaN;
ratemap_down_left_rotated = ratemap_comp2;

%ratemp comp 3 rotated 180°
[posrx_comp3,posry_comp3] = center_path(posrx(bottomright),posry(bottomright),shape);
posrx_comp3=-posrx_comp3;posry_comp3=-posry_comp3;
[mapAxis_comp3] = mapaxis(posrx_comp3,posry_comp3,binWidth);
visited_comp3 = visitedBins(posrx_comp3,posry_comp3,mapAxis_comp3); % calculate the visited zones of the arena 
[spkx_comp3,spky_comp3] = get_pos_spikes(ts(spike_bottom_right),posrx_comp3,posry_comp3,posts(bottomright)); % get the position of the spikes
[ratemap_comp3] = ratemap_gaussian(2*binWidth,spkx_comp3,spky_comp3,posrx_comp3,posry_comp3,posts(bottomright),binWidth,mapAxis_comp3); % build matrix ratemap 
ratemap_comp3(visited_comp3==0) = NaN;
ratemap_down_right_rotated = ratemap_comp3;

%ratemp comp 4 rotated -90°
[posrx_comp4,posry_comp4] = center_path(posry(topright),posrx(topright),shape);
posrx_comp4=-posrx_comp4;
[mapAxis_comp4] = mapaxis(posrx_comp4,posry_comp4,binWidth);
visited_comp4 = visitedBins(posrx_comp4,posry_comp4,mapAxis_comp4); % calculate the visited zones of the arena 
[spkx_comp4,spky_comp4] = get_pos_spikes(ts(spike_top_right),posrx_comp4,posry_comp4,posts(topright)); % get the position of the spikes
[ratemap_comp4] = ratemap_gaussian(2*binWidth,spkx_comp4,spky_comp4,posrx_comp4,posry_comp4,posts(topright),binWidth,mapAxis_comp4); % build matrix ratemap 
ratemap_comp4(visited_comp4==0) = NaN;
ratemap_top_right_rotated = ratemap_comp4;

minLength = min([length(ratemap_top_left2), length(ratemap_down_left_rotated), length(ratemap_down_right_rotated), length(ratemap_top_right_rotated)]);
ratemap_top_left2=ratemap_top_left2(1:minLength,1:minLength); ratemap_down_left_rotated=ratemap_down_left_rotated(1:minLength,1:minLength); ratemap_down_right_rotated=ratemap_down_right_rotated(1:minLength,1:minLength); ratemap_top_right_rotated=ratemap_top_right_rotated(1:minLength,1:minLength); 

[corrcoeff1_rotated,nonZero_corrcoeff1_rotated,RrankO1_rotated] = get_2D_correlation_value(ratemap_top_left2,ratemap_down_left_rotated);
[corrcoeff2_rotated,nonZero_corrcoeff2_rotated,RrankO2_rotated] = get_2D_correlation_value(ratemap_top_left2,ratemap_down_right_rotated);
[corrcoeff3_rotated,nonZero_corrcoeff3_rotated,RrankO3_rotated] = get_2D_correlation_value(ratemap_top_left2,ratemap_top_right_rotated);
[corrcoeff4_rotated,nonZero_corrcoeff4_rotated,RrankO4_rotated] = get_2D_correlation_value(ratemap_down_left_rotated,ratemap_down_right_rotated);
[corrcoeff5_rotated,nonZero_corrcoeff5_rotated,RrankO5_rotated] = get_2D_correlation_value(ratemap_down_left_rotated,ratemap_top_right_rotated);
[corrcoeff6_rotated,nonZero_corrcoeff6_rotated,RrankO6_rotated] = get_2D_correlation_value(ratemap_down_right_rotated,ratemap_top_right_rotated);

meancorrcoef_rotated = mean([corrcoeff1_rotated,corrcoeff2_rotated,corrcoeff3_rotated,corrcoeff4_rotated,corrcoeff5_rotated,corrcoeff6_rotated]);

four_comp_corr = figure('visible', 'off');
subplot(2,5,1);
myfig=pcolor(ratemap_top_left2);colormap jet;set(myfig,'Edgecolor','none');axis xy;set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);axis off;
title(strcat('mean spatial corr = ',num2str(meancorrcoef)));
subplot(2,5,2);
myfig=pcolor(ratemap_top_right2);colormap jet;set(myfig,'Edgecolor','none');axis xy;set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);axis off;
subplot(2,5,6);
myfig=pcolor(ratemap_down_left2);colormap jet;set(myfig,'Edgecolor','none');axis xy;set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);axis off;
subplot(2,5,7);
myfig=pcolor(ratemap_down_right2);colormap jet;set(myfig,'Edgecolor','none');axis xy;set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);axis off;
subplot(2,5,4);
myfig=pcolor(ratemap_top_left2);colormap jet;set(myfig,'Edgecolor','none');axis xy;set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);axis off;
title(strcat('mean spatial corr = ',num2str(meancorrcoef_rotated)));
subplot(2,5,5);
myfig=pcolor(ratemap_top_right_rotated);colormap jet;set(myfig,'Edgecolor','none');axis xy;set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);axis off;
subplot(2,5,9);
myfig=pcolor(ratemap_down_left_rotated);colormap jet;set(myfig,'Edgecolor','none');axis xy;set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);axis off;
subplot(2,5,10);
myfig=pcolor(ratemap_down_right_rotated);colormap jet;set(myfig,'Edgecolor','none');axis xy;set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);axis off;

set(four_comp_corr,'units','normalized','outerposition',[0 0 1 1 ]);
set(four_comp_corr,'Units','Inches');
pos = get(four_comp_corr,'Position');
set(four_comp_corr,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

   if save_figure == 1
            complete_name = strcat(pathfigure,inputfile(1:end-4),'_t',num2str(p),'_c',num2str(i),'_spatial_correlation_four_comp');
            print('-dpng',complete_name,'-r100');
%             print(two_comp_corr,complete_name,'-dpdf')

   end
   if close_figure == 1
        close(four_comp_corr);
   end

%%
posts_topleft = posts(topleft);
posts_down_left = posts(bottomleft);
posts_top_right = posts(topright);
posts_down_right = posts(bottomright);

diffposts_topleft = diff(posts_topleft);
diffsuperioto500ms = find(diffposts_topleft>500);
rate_topleft = [];duration_topleft = [];
for hgt = 1:length(diffsuperioto500ms)-1
%     path_x = posrx_in(diffsuperioto500ms(hgt)+1:diffsuperioto500ms(hgt+1));
%     path_y = posry_in(diffsuperioto500ms(hgt)+1:diffsuperioto500ms(hgt+1));
    time = posts_topleft(diffsuperioto500ms(hgt)+1:diffsuperioto500ms(hgt+1));
    spike = find(ts>posts_topleft(diffsuperioto500ms(hgt)+1) & ts<posts_topleft(diffsuperioto500ms(hgt+1)));
    duration = (time(end)-time(1))/1000;
    numspike = length(spike);
    rate_topleft(hgt)=numspike/duration;
    duration_topleft(hgt) = duration;

end
rate_topleft_without0= rate_topleft;
rate_topleft_without0(rate_topleft_without0==0) = [];

diffposts_down_left = diff(posts_down_left);
diffsuperioto500ms = find(diffposts_down_left>500);
rate_down_left = [];duration_down_left = [];
for hgt = 1:length(diffsuperioto500ms)-1
%     path_x = posrx_in(diffsuperioto500ms(hgt)+1:diffsuperioto500ms(hgt+1));
%     path_y = posry_in(diffsuperioto500ms(hgt)+1:diffsuperioto500ms(hgt+1));
    time = posts_down_left(diffsuperioto500ms(hgt)+1:diffsuperioto500ms(hgt+1));
    spike = find(ts>posts_down_left(diffsuperioto500ms(hgt)+1) & ts<posts_down_left(diffsuperioto500ms(hgt+1)));
    duration = (time(end)-time(1))/1000;
    numspike = length(spike);
    rate_down_left(hgt)=numspike/duration;
    duration_down_left(hgt) = duration;

end
rate_down_left_without0= rate_down_left;
rate_down_left_without0(rate_down_left_without0==0) = [];

diffposts_top_right = diff(posts_top_right);
diffsuperioto500ms = find(diffposts_top_right>500);
rate_top_right = [];duration_top_right = [];
for hgt = 1:length(diffsuperioto500ms)-1
%     path_x = posrx_in(diffsuperioto500ms(hgt)+1:diffsuperioto500ms(hgt+1));
%     path_y = posry_in(diffsuperioto500ms(hgt)+1:diffsuperioto500ms(hgt+1));
    time = posts_top_right(diffsuperioto500ms(hgt)+1:diffsuperioto500ms(hgt+1));
    spike = find(ts>posts_top_right(diffsuperioto500ms(hgt)+1) & ts<posts_top_right(diffsuperioto500ms(hgt+1)));
    duration = (time(end)-time(1))/1000;
    numspike = length(spike);
    rate_top_right(hgt)=numspike/duration;
    duration_top_right(hgt) = duration;

end
rate_top_right_without0= rate_top_right;
rate_top_right_without0(rate_top_right_without0==0) = [];


diffposts_down_right = diff(posts_down_right);
diffsuperioto500ms = find(diffposts_down_right>500);
rate_down_right = [];duration_down_right = [];
for hgt = 1:length(diffsuperioto500ms)-1
%     path_x = posrx_in(diffsuperioto500ms(hgt)+1:diffsuperioto500ms(hgt+1));
%     path_y = posry_in(diffsuperioto500ms(hgt)+1:diffsuperioto500ms(hgt+1));
    time = posts_down_right(diffsuperioto500ms(hgt)+1:diffsuperioto500ms(hgt+1));
    spike = find(ts>posts_down_right(diffsuperioto500ms(hgt)+1) & ts<posts_down_right(diffsuperioto500ms(hgt+1)));
    duration = (time(end)-time(1))/1000;
    numspike = length(spike);
    rate_down_right(hgt)=numspike/duration;
    duration_down_right(hgt) = duration;

end
rate_down_right_without0= rate_down_right;
rate_down_right_without0(rate_down_right_without0==0) = [];

%% plots for 4 compartments and calculate path activity   

temporal_four_comp = figure('visible', 'off');
spkx_topleft = spkx(spike_top_left);
spky_topleft = spky(spike_top_left);
subplot(3,4,1);
plot(posx(topleft),posy(topleft),'color',[.4 .4 .4],'linewidth',2);%,spkx,spky,'.r');
 hold on
scatter(spkx_topleft,spky_topleft,400,'.r');
set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);
axis off
box off

spkx_downleft = spkx(spike_bottom_left);
spky_downleft = spky(spike_bottom_left);
subplot(3,4,2);
plot(posx(bottomleft),posy(bottomleft),'color',[.4 .4 .4],'linewidth',2);%,spkx,spky,'.r');
hold on
scatter(spkx_downleft,spky_downleft,400,'.r');
set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);
axis off
box off

spkx_topright = spkx(spike_top_right);
spky_topright = spky(spike_top_right);
subplot(3,4,3);
plot(posx(topright),posy(topright),'color',[.4 .4 .4],'linewidth',2);%,spkx,spky,'.r');
 hold on
scatter(spkx_topright,spky_topright,400,'.r');
set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);
axis off
box off

spkx_downright = spkx(spike_bottom_right);
spky_downright = spky(spike_bottom_right);
subplot(3,4,4);
plot(posx(bottomright),posy(bottomright),'color',[.4 .4 .4],'linewidth',2);%,spkx,spky,'.r');
 hold on
scatter(spkx_downright,spky_downright,400,'.r');
set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);
axis off
box off

subplot(3,4,5);
plot(rate_topleft,'k-o','LineWidth',1,'MarkerSize',5,'MarkerEdgeColor','k','MarkerFaceColor',[0.5,0.5,0.5]);
axis square, box off;
h=get(gcf, 'currentaxes');
set(h, 'fontsize', 10, 'linewidth', 2);
box off;
get(gca,'ticklength');
set(gca,'TickLength',[0.03,0.02]);
set(gca,'TickDir','out');
ylim([0 max([rate_topleft,rate_down_left,rate_top_right,rate_down_right])]);
xlabel('path'); ylabel('frequency(Hz)');

subplot(3,4,6);
plot(rate_down_left,'k-o','LineWidth',1,'MarkerSize',5,'MarkerEdgeColor','k','MarkerFaceColor',[0.5,0.5,0.5]);
axis square, box off;
h=get(gcf, 'currentaxes');
set(h, 'fontsize', 10, 'linewidth', 2);
box off;
get(gca,'ticklength');
set(gca,'TickLength',[0.03,0.02]);
set(gca,'TickDir','out');
ylim([0 max([rate_topleft,rate_down_left,rate_top_right,rate_down_right])]);
xlabel('path'); ylabel('frequency(Hz)');

subplot(3,4,7);
plot(rate_top_right,'k-o','LineWidth',1,'MarkerSize',5,'MarkerEdgeColor','k','MarkerFaceColor',[0.5,0.5,0.5]);
axis square, box off;
h=get(gcf, 'currentaxes');
set(h, 'fontsize', 10, 'linewidth', 2);
box off;
get(gca,'ticklength');
set(gca,'TickLength',[0.03,0.02]);
set(gca,'TickDir','out');
ylim([0 max([rate_topleft,rate_down_left,rate_top_right,rate_down_right])]);
xlabel('path'); ylabel('frequency(Hz)');

subplot(3,4,8);
plot(rate_down_right,'k-o','LineWidth',1,'MarkerSize',5,'MarkerEdgeColor','k','MarkerFaceColor',[0.5,0.5,0.5]);
axis square, box off;
h=get(gcf, 'currentaxes');
set(h, 'fontsize', 10, 'linewidth', 2);
box off;
get(gca,'ticklength');
set(gca,'TickLength',[0.03,0.02]);
set(gca,'TickDir','out');
ylim([0 max([rate_topleft,rate_down_left,rate_top_right,rate_down_right])]);
xlabel('path'); ylabel('frequency(Hz)');

subplot(3,4,9);
plot(rate_topleft_without0,'k-o','LineWidth',2,'MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor',[0.5,0.5,0.5]);
axis square, box off;
h=get(gcf, 'currentaxes');
set(h, 'fontsize', 15, 'linewidth', 2);
box off;
get(gca,'ticklength');
set(gca,'TickLength',[0.03,0.02]);
set(gca,'TickDir','out');
ylim([0 max([rate_topleft_without0,rate_down_left_without0,rate_top_right_without0,rate_down_right_without0])]);
xlabel('path'); ylabel('frequency(Hz)');

subplot(3,4,10);
plot(rate_down_left_without0,'k-o','LineWidth',2,'MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor',[0.5,0.5,0.5]);
axis square, box off;
h=get(gcf, 'currentaxes');
set(h, 'fontsize', 15, 'linewidth', 2);
box off;
get(gca,'ticklength');
set(gca,'TickLength',[0.03,0.02]);
set(gca,'TickDir','out');
ylim([0 max([rate_topleft_without0,rate_down_left_without0,rate_top_right_without0,rate_down_right_without0])]);
xlabel('path'); ylabel('frequency(Hz)');

subplot(3,4,11);
plot(rate_top_right_without0,'k-o','LineWidth',2,'MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor',[0.5,0.5,0.5]);
axis square, box off;
h=get(gcf, 'currentaxes');
set(h, 'fontsize', 15, 'linewidth', 2);
box off;
get(gca,'ticklength');
set(gca,'TickLength',[0.03,0.02]);
set(gca,'TickDir','out');
ylim([0 max([rate_topleft_without0,rate_down_left_without0,rate_top_right_without0,rate_down_right_without0])]);
xlabel('path'); ylabel('frequency(Hz)');

subplot(3,4,12);
plot(rate_down_right_without0,'k-o','LineWidth',2,'MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor',[0.5,0.5,0.5]);
axis square, box off;
h=get(gcf, 'currentaxes');
set(h, 'fontsize', 15, 'linewidth', 2);
box off;
get(gca,'ticklength');
set(gca,'TickLength',[0.03,0.02]);
set(gca,'TickDir','out');
ylim([0 max([rate_topleft_without0,rate_down_left_without0,rate_top_right_without0,rate_down_right_without0])]);
xlabel('path'); ylabel('frequency(Hz)');

set(temporal_four_comp,'units','normalized','outerposition',[0 0 1 1 ]);
set(temporal_four_comp,'Units','Inches');
pos = get(temporal_four_comp,'Position');
set(temporal_four_comp,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])


   if save_figure == 1
            complete_name = strcat(pathfigure,inputfile(1:end-4),'_t',num2str(p),'_c',num2str(i),'_temporal_activity_four_comp');
            print('-dpng',complete_name,'-r100');

   end
   if close_figure == 1
        close(temporal_four_comp);
   end

%% temporal activity with all compartments next to the other
temporal_four_comp_all = figure('visible', 'off');
path_order = zeros(length(posts),1);
path_order(topleft)=10;
path_order(topright)=-10;
path_order(bottomright)=20;
path_order(bottomleft)=-20;

path_order(door_top)=1;
path_order(door_bottom)=1;
path_order(door_left)=1;
path_order(door_right)=1;
path_order(path_order==0)=1000;

if path_order(end)==1000
        gh = length(path_order);
        while path_order(gh) == 1000
            path_order(gh) = [];
            gh = gh-1;
        end    
end

copy_path = path_order;
for ii = 1:length(path_order)

    if path_order(ii) == 1000
        j = ii;
        while path_order(j) == 1000
            j = j+1;
        end

        next_val = path_order(j);
        copy_path(ii) = next_val;


    end
end

transition = abs(diff(copy_path));
[id,~] = find(transition>0);

remove_short_id = diff(id);
id(remove_short_id==1)=[];

rate = [];duration_path = [];

    time = posts(1:id(1));
    spike = find(ts>posts(1) & ts<posts(id(1)));
    duration = (time(end)-time(1))/1000;
    numspike = length(spike);
    rate(1,1)=numspike/duration;
    duration_path(1,1) = duration;

for thtg = 1:length(id)-1
    time = posts(id(thtg):id(thtg+1));
    spike = find(ts>posts(id(thtg)) & ts<posts(id(thtg+1)));
    duration = (time(end)-time(1))/1000;
    numspike = length(spike);
    rate(thtg+1)=numspike/duration;
    duration_path(thtg+1) = duration;   
end

    time = posts(id(end):end);
    spike = find(ts>posts(id(end)) & ts<posts(end));
    duration = (time(end)-time(1))/1000;
    numspike = length(spike);
    rate(1,end+1)=numspike/duration;
    duration_path(1,end+1) = duration;

compartment_side = copy_path(id);
rate = rate';
rate(rate == 0) = -1;
idx_top_left = find(compartment_side == 10);
idx_top_right = find(compartment_side == -10);
idx_bottom_left = find(compartment_side == -20);
idx_bottom_right = find(compartment_side == 20); 

scatter(idx_top_left,rate(compartment_side == 10),150,'or','filled');
hold on
scatter(idx_top_right,rate(compartment_side == -10),150,'ob','filled');
scatter(idx_bottom_left,rate(compartment_side == -20),150,'og','filled');
scatter(idx_bottom_right,rate(compartment_side == 20),150,'o','filled','MarkerEdgeColor',[1 0.75 0],'MarkerfaceColor',[1 0.75 0]);
box off;
h=get(gcf, 'currentaxes');
set(h, 'fontsize', 15, 'linewidth', 2);
box off;
get(gca,'ticklength');
set(gca,'TickLength',[0.03,0.02]);
set(gca,'TickDir','out');
ylim([0 max(rate)]);
xlabel('path'); ylabel('frequency(Hz)');
legend({'top left','top right','bottom left','bottom right'})%,'Position',[.1 .1 .2 .2]);

set(temporal_four_comp_all,'units','normalized','outerposition',[0 0 1 1 ]);
set(temporal_four_comp_all,'Units','Inches');
pos = get(temporal_four_comp_all,'Position');
set(temporal_four_comp_all,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])


   if save_figure == 1
            complete_name = strcat(pathfigure,inputfile(1:end-4),'_t',num2str(p),'_c',num2str(i),'_temporal_activity_all_four_comp');
            print('-dpng',complete_name,'-r100');

   end
   if close_figure == 1
        close(temporal_four_comp_all);
   end


%% save results
results_txt_4comp(results_folder,inputfile,p,i,...
    peak_rate_global,peak_rate_top_left,peak_rate_bottom_left,peak_rate_bottom_right,peak_rate_top_right,...
    peak_rate_door_top,peak_rate_door_left,peak_rate_door_bottom,peak_rate_door_right,...
    meancorrcoef,meancorrcoef_rotated,corrcoeff1,corrcoeff2,corrcoeff3,corrcoeff4,corrcoeff5,corrcoeff6,...
    corrcoeff1_rotated,corrcoeff2_rotated,corrcoeff3_rotated,corrcoeff4_rotated,corrcoeff5_rotated,corrcoeff6_rotated,...
    ret,peak_frequency_theta,theta_power,delta_power,sig_threshold,theta_modulation,...
    PFD_deg,Vector_length,peak_HD,...
    PFD_deg_top_left,Vector_length_top_left,peak_HD_top_left,...
    PFD_deg_bottom_left,Vector_length_bottom_left,peak_HD_bottom_left,...
    PFD_deg_bottom_right,Vector_length_bottom_right,peak_HD_bottom_right,...
    PFD_deg_top_right,Vector_length_top_right,peak_HD_top_right,...
    PFD_deg_door_top,Vector_length_door_top,peak_HD_door_top,...
    PFD_deg_door_left,Vector_length_door_left,peak_HD_door_left,...
    PFD_deg_door_bottom,Vector_length_door_bottom,peak_HD_door_bottom,...
    PFD_deg_door_right,Vector_length_door_right,peak_HD_door_right,...
    Vector_length_double,flipping_score,flipping_score_range,Vector_length_quadri,Four_dir_score,Four_dir_score_minmax,...
    Four_dir_score_range,Four_dir_score_minmax_range,number_of_corr_peaks);

%% save mat
nameformat = strcat(mat_folder,filename(1:end-4),'_t',num2str(p),'_c',num2str(i),'.mat');
save(nameformat,...
    'TTMtx2','PosMtx','filename','ts','waveforms','xy_top_left','xy_bottom_left','xy_top_right','xy_bottom_right',...
    'xy_top_door','xy_left_door','xy_bottom_door','xy_right_door',...
    'topleft','bottomleft','bottomright','topright',...
    'topdoor','leftdoor','bottomdoor','rightdoor',...
    'spike_top_left','spike_bottom_left','spike_top_right','spike_bottom_right',...
    'door_top','door_left','door_bottom','door_right',...
    'ratemap_top_left','ratemap_bottom_left', 'ratemap_bottom_right','ratemap_top_right',...
    'ratemap_door_top','ratemap_door_left', 'ratemap_door_bottom','ratemap_door_right',...
    'ratemap_top_left2','ratemap_down_left2','ratemap_down_right2','ratemap_top_right2',...
    'ratemap_down_left_rotated','ratemap_down_right_rotated','ratemap_top_right_rotated',...
    'dirrate','dirrate_top_left','dirrate_bottom_left','dirrate_bottom_right','dirrate_top_right',...
    'dirrate_door_top','dirrate_door_left','dirrate_door_bottom','dirrate_door_right', ...
    'corrcoeff1','corrcoeff2','corrcoeff3','corrcoeff4','corrcoeff5','corrcoeff6',...
    'corrcoeff1_rotated','corrcoeff2_rotated','corrcoeff3_rotated','corrcoeff4_rotated','corrcoeff5_rotated','corrcoeff6_rotated');
%%

end
    end
    end
%     end
    
        end
end
