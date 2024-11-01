%% Figure peak difference between comp

cd 'C:\Users\celia\Pictures\ConnecSpace\FINAL ANALYSIS\Result\Protocole flip cell\';
[num,txt,raw] = xlsread('Suivi directional cells','BDC properties');

% S1
% Inputfile = string(txt(2:28,2));
% Tetrode = string(raw(2:28,3));
% Cell = string(raw(2:28,4));
% session = "S1";

% S2
% Inputfile = string(txt(35:69,2));
% Tetrode = string(raw(35:69,3));
% Cell = string(raw(35:69,4));
% session = "S2";

% S3
Inputfile = string(txt(94:125,2));
Tetrode = string(raw(94:125,3));
Cell = string(raw(94:125,4));
session = "S3";

mat_file = 'C:\Users\celia\Pictures\ConnecSpace\FINAL ANALYSIS\Mat\Protocoles flip cell\S3\';
mat_file2 = 'C:\Users\celia\Pictures\ConnecSpace\FINAL ANALYSIS\Mat\Protocoles flip cell\Protocole flip cell - S1 et S2 inversés\S3 4 comp\';
mat_file3 = 'C:\Users\celia\Pictures\ConnecSpace\FINAL ANALYSIS\Mat\Protocoles flip cell\Cue même côté S1 et S2\S3\';

% define paramaters:
binSizeDir = 6;    % bining for polar plot (choose 5,6 or 10)
pathfigure = 'C:\Users\celia\Pictures\ConnecSpace\FINAL ANALYSIS\Result\Protocole flip cell\Peak between compartment\'; % folder to save figure (finish the path with\)
diff_load = 0;

for i = 1:size(Inputfile,1)
    cd(mat_file);
        
    filename = Inputfile{i};
    diff_inputfile = size(strcat(filename(1:end-4),'_t',Tetrode{i},'_c',Cell{i},'.mat'),2)-size(filename,2);
    
    inputfile(i,1:size(filename,2)+diff_inputfile) = strcat(filename(1:end-4),'_t',Tetrode{i},'_c',Cell{i},'.mat');
    
    if i > 1
        diff_load = abs(size(inputfile,2)-size(strcat(filename(1:end-4),'_t',Tetrode{i},'_c',Cell{i},'.mat'),2));   
    end
    
    if exist(inputfile(i,:)) == 2
        load(inputfile(i,1:end-diff_load));
    else 
        cd(mat_file2);
        if exist(inputfile(i,:)) == 2
            load(inputfile(i,1:end-diff_load));
        else
            cd(mat_file3);
            load(inputfile(i,1:end-diff_load));        
        end
    end
    % dirrate_top_right(60)=mean([dirrate_top_right(1),dirrate_top_right(59)]);
    
    % title_file=replace(inputfile(i,1:25),'_','-');
    fig = figure('units','normalized','outerposition',[0 0 1 1]);
    sgtitle({
      session
    %   title_file
      },'FontSize',30);
    
    if session == "S3"
        subplot(3,2,1:2)
    else
        subplot(2,2,1:2)
    end
    
    line([0 0],[-max(dirrate) max(dirrate)],'linewidth',2,'color','k');
    hold on
    line([-max(dirrate) max(dirrate)],[0 0],'linewidth',2,'color','k');
    Angles = (0:binSizeDir:360-binSizeDir)';
    nbins = numel(Angles);
    x = cosd(Angles).*dirrate;
    y = sind(Angles).*dirrate;
    X = smooth(x,3);
    Y = smooth(y,3);
    hold on ;
    v = [X Y];
    f = (1:1:length(dirrate));
    cc = colormapc(7,nbins);
    col = cc;
    patch('Faces',f,'Vertices',v,'FaceVertexCData',col,'EdgeColor','interp','FaceColor','non','LineWidth',3,'LineSmoothing','off');
    hold off;
    axis equal tight;
    axis off;
    
    peak_HD = max(dirrate);
    peak_HD_2=sprintf('%.2f', peak_HD);
    set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);
    subtitle({['Global peak = ', peak_HD_2,'Hz']},'fontsize',20);
    set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);xlabel('global')
    
    if session == "S3"
        subplot(3,2,3)
        dirrate_left=dirrate_top_left;
        dirrate_right=dirrate_top_right;
    else
        subplot(2,2,3)
    end
    
    Rate_max = max(max(dirrate_left),max(dirrate_right));
    if Rate_max == max(dirrate_left)
        MaxBinRate = find(dirrate_left == max(dirrate_left)) ; 
    else
        MaxBinRate = find(dirrate_right == max(dirrate_right)) ; 
    end
    
    line([0 0],[-(Rate_max) Rate_max],'linewidth',2,'color','k');
    hold on
    line([-Rate_max Rate_max],[0 0],'linewidth',2,'color','k');
    cc = colormapc(5,nbins);
    ColorForMaxPeak = [1 0 0]   ; % This is Blue but you can change to whatever colour you specify from colormapc function
    MaxBinColorMap = find(ismember(cc,ColorForMaxPeak,'rows'))   ;
    MaxBinColorMap=MaxBinColorMap(1);
    ToShift = MaxBinRate(1) -  MaxBinColorMap ;
    cc = circshift(cc,ToShift);
    
    x = cosd(Angles).*dirrate_left;
    y = sind(Angles).*dirrate_left;
    X = smooth(x,3);
    Y = smooth(y,3);
    hold on ;
    v = [X Y];
    f = (1:1:length(dirrate_left));
    col = cc;
    patch('Faces',f,'Vertices',v,'FaceVertexCData',col,'EdgeColor','interp','FaceColor','non','LineWidth',3,'LineSmoothing','off');
    hold off;
    axis equal tight;
    axis off;
    
    peak_HD = max(dirrate_left);
    peak_HD_2=sprintf('%.2f', peak_HD);
    set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);
    subtitle({['Left peak = ', peak_HD_2,'Hz']},'fontsize',20);
    set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);xlabel('global')
    
    if session == "S3"
        subplot(3,2,4)
    else
        subplot(2,2,4)
    end
    
    line([0 0],[-(Rate_max) Rate_max],'linewidth',2,'color','k');
    hold on
    line([-Rate_max Rate_max],[0 0],'linewidth',2,'color','k');
    Angles = (0:binSizeDir:360-binSizeDir)';
    nbins = numel(Angles);
    
    x = cosd(Angles).*dirrate_right;
    y = sind(Angles).*dirrate_right;
    X = smooth(x,3);
    Y = smooth(y,3);
    hold on ;
    v = [X Y];
    f = (1:1:length(dirrate_right));
    col = cc;
    patch('Faces',f,'Vertices',v,'FaceVertexCData',col,'EdgeColor','interp','FaceColor','non','LineWidth',3,'LineSmoothing','off');
    hold off;
    axis equal tight;
    axis off;
    
    peak_HD = max(dirrate_right);
    peak_HD_2=sprintf('%.2f', peak_HD);
    set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);
    subtitle({['Right peak = ', peak_HD_2,'Hz']},'fontsize',20);
    set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);xlabel('global')
    
    if session == "S3"
        subplot(3,2,5)
    
        line([0 0],[-(Rate_max) Rate_max],'linewidth',2,'color','k');
        hold on
        line([-Rate_max Rate_max],[0 0],'linewidth',2,'color','k');
        Angles = (0:binSizeDir:360-binSizeDir)';
        nbins = numel(Angles);
        
        x = cosd(Angles).*dirrate_bottom_left;
        y = sind(Angles).*dirrate_bottom_left;
        X = smooth(x,3);
        Y = smooth(y,3);
        hold on ;
        v = [X Y];
        f = (1:1:length(dirrate_bottom_left));
        col = cc;
        patch('Faces',f,'Vertices',v,'FaceVertexCData',col,'EdgeColor','interp','FaceColor','non','LineWidth',3,'LineSmoothing','off');
        % circle(0,0,Rate_max);
        hold off;
        axis equal tight;
        axis off;
        
        peak_HD = max(dirrate_bottom_left);
        peak_HD_2=sprintf('%.2f', peak_HD);
        set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);
        subtitle({['Right peak = ', peak_HD_2,'Hz']},'fontsize',20);
        set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);xlabel('global')
        
        
        subplot(3,2,6)
        
        
        line([0 0],[-(Rate_max) Rate_max],'linewidth',2,'color','k');
        hold on
        line([-Rate_max Rate_max],[0 0],'linewidth',2,'color','k');
        Angles = (0:binSizeDir:360-binSizeDir)';
        nbins = numel(Angles);
        
        x = cosd(Angles).*dirrate_bottom_right;
        y = sind(Angles).*dirrate_bottom_right;
        X = smooth(x,3);
        Y = smooth(y,3);
        hold on ;
        v = [X Y];
        f = (1:1:length(dirrate_bottom_right));
        col = cc;
        patch('Faces',f,'Vertices',v,'FaceVertexCData',col,'EdgeColor','interp','FaceColor','non','LineWidth',3,'LineSmoothing','off');
        % circle(0,0,Rate_max);
        hold off;
        axis equal tight;
        axis off;
        
        peak_HD = max(dirrate_bottom_right);
        peak_HD_2=sprintf('%.2f', peak_HD);
        set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);
        subtitle({['Right peak = ', peak_HD_2,'Hz']},'fontsize',20);
        set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);xlabel('global')
    end

    complete_name = strcat(pathfigure,inputfile(i,1:end-diff_load-4));
%     print('-dpng',complete_name,'-r500');
    
    close();

end