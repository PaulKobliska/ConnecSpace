% *************************************************************************
% script to analyze 2 and 4 compartments HD activities
% handmade script 2019 M2 internship Célia
% corrections made by P-Yves 15/03/2019
%modification made by Celia 05/01/2021
%modification made by Célia 05/09/2022 : change selection of compartments and
%door by crating polygon with mouse clics.
% *************************************************************************

pathfigure = ''; % enter folder to save figure (finish the path with\)
results_folder = ''; % enter folder to save data

a= dir; 
b= struct2table(a);
all_files= b.name(:);
% remove directories from files list
all_files([1,2],:)=[];
zz = size (all_files,1);

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

%%
for f = 1:zz

    filename = all_files(f);
    filename =filename{1};
    inputfile = filename;
    
    % *************************************************************************
    shape = 'circle';
    
    load(inputfile);
    
    [posx,posy,posts,mapAxis,visited] = posdata_multicompartment(PosMtx,threshold,shape,binWidth);
    
    [xy_left,xy_right,xy_door,left,right,door] = Splitting_compartment_v1_removing_door(PosMtx,inputfile);
    close
    
    TTMtx=TTMtx2;
    tetrode = unique(TTMtx(:,2));    
    for pp = 1:length(tetrode)
        TTMtx2=TTMtx(TTMtx(:,2)== tetrode(pp),:);
        TTMtx2(TTMtx2(:,3)==0,:) = [];
        p = tetrode(pp);
        cell = unique(TTMtx2(:,3));
        for iii = 1:length(cell)
            ts = TTMtx2(TTMtx2(:,3)==cell(iii),1);
            cells = cell(iii);
            waveforms = TTMtx2(TTMtx2(:,3)==cell(iii),4:131);
            i = cells;
    
            %%
            [spkx,spky,spike_left,spike_right,spike_door] = Splitting_two_compartment_and_door(ts,posx,posy,posts,xy_left,xy_right,xy_door);
        
            two_comp = figure('visible', 'off');
         
            %% Rate map global
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
           end
        
            if close_figure == 1
                close(two_comp);
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
            [rotated_mapAxis_comp1] = mapAxis_comp1;
            rotated_visited_comp1 = visitedBins(minus_posrx_comp1,minus_posry_comp1,rotated_mapAxis_comp1); % calculate the visited zones of the arena 
            [rotated_spkx_comp1,rotated_spky_comp1] = get_pos_spikes(ts(spike_left),minus_posrx_comp1,minus_posry_comp1,posts(left)); % get the position of the spikes
            [rotated_ratemap_comp1] = ratemap_gaussian(2*binWidth,rotated_spkx_comp1,rotated_spky_comp1,minus_posrx_comp1,minus_posry_comp1,posts(left),binWidth,rotated_mapAxis_comp1); % build matrix ratemap 
            rotated_ratemap_comp1(rotated_visited_comp1==0) = NaN;
            
            % ratemap comp 2
            [posrx_comp2,posry_comp2] = center_path(posrx(right),posry(right),shape);
            [mapAxis_comp2] = mapAxis_comp1;
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
           end
           if close_figure == 1
                close(colored_plot);
           end
            
            %% Save results
            results_txt_2comp(results_folder,inputfile,p,i,peak_rate_global,peak_rate_left,peak_rate_right,...
                peak_rate_door,corrcoeff,rotated_corrcoeff,PFD_deg,Vector_length,peak_HD,PFD_deg_left,...
                Vector_length_left,peak_HD_left,PFD_deg_right,Vector_length_right,peak_HD_right,...
                PFD_deg_door,Vector_length_door,peak_HD_door,Vector_length_double,flipping_score,...
                flipping_score_range,Vector_length_quadri,Four_dir_score,Four_dir_score_minmax,...
                Four_dir_score_range,Four_dir_score_minmax_range,number_of_corr_peaks);
        
        end
    end
end
