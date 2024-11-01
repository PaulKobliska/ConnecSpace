% function  polarplot_v2(OpenfIeld,posgx,posgy,posrx,posry,posts,ts,threshold,shape)
function [PFD_deg,Vector_length,peak_HD,dirrate,deg_theta,posdirts,dirmap]= polarplot_compartment_new(posts,posgx,posgy,posrx,posry,ts,correction,binSizeDir)
%polar coordinates
samplerate = 25;

thetacartesian = mod((180/pi)*(atan2(-posry+posgy, posrx-posgx)),360);
% thetacartesian2 = (atan2(-posry+posgy, posrx-posgx));
deg_theta = round(thetacartesian);
% deg_theta = smooth(deg_theta,5);

posdirts = zeros(length(posts),1);
posdirts(1,1) = sum(ts > 0 & ts <= posts(1,1));
for l = 2:size(posdirts,1)
    posdirts(l,1) = sum(ts > posts(l-1,1) & ts <= posts(l,1));
end

% PosMtx(:,6) = 0;
% PosMtx(1,6) = sum(ts > 0 & ts <= PosMtx(1,1));
% for l = 2:size(PosMtx,1)
%     PosMtx(l,6) = sum(ts > PosMtx(l-1,1) & ts <= PosMtx(l,1));
% end

degbin=0:binSizeDir:360;
dirmap = zeros(size(degbin,2)-1,2);
for b = 1:size(degbin,2)-1
    dirmap(b,1) = sum(deg_theta > degbin(b) & deg_theta <= degbin(b+1));
    dirmap(b,2) = sum(posdirts(deg_theta > degbin(b) & deg_theta <= degbin(b+1),1));
end
% %corrections
% factor = length(dirmap)/4;
% dirmap_v2 = repmat(dirmap,3);
% 
% dirmap(16,:)= round(mean(dirmap(15:17,:)));
% dirmap(31,:)= round(mean(dirmap(30:32,:)));
% dirmap(46,:)= round(mean(dirmap(45:47,:)));
% dirmap(1,:)= round(mean(dirmap(1:2,:)));
% dirmap(60,:)= round(mean(dirmap(59:60,:)));

dirrate=dirmap(:,2)./(dirmap(:,1)/samplerate);
dirrate(isnan(dirrate))=0;
% dirrate(end,1) =dirrate(1,1);
dirrate = smooth(dirrate,3);
%correction of animal head relative to LEDs

dirrate = circshift(dirrate,-correction/binSizeDir);

Angles = [0:binSizeDir:360-binSizeDir]';
nbins = numel(Angles);
Rate = [dirrate];
x = cosd(Angles).*Rate;
y = sind(Angles).*Rate;

% MaxBinRate = find(Rate == max(Rate)) ; 
% MaxAngleFiringRate = Angles(MaxBinRate) ; 
cc = colormapc(7,nbins);
%h = figure(3);
% polar plot in blue
%ColorForMaxPeak = [0 0 1]   ; % This is Blue but you can change to whatever colour you specify from colormapc function
% ColorForMaxPeak = [1 0 0] ; %red color
% ColorForMaxPeak = [0 1 0] ; %green color
% ColorForMaxPeak = [0 1 0];

% MaxBinColorMap = find(ismember(cc,ColorForMaxPeak,'rows'))   ;MaxBinColorMap=MaxBinColorMap(1);
% ToShift = MaxBinRate -  MaxBinColorMap ;
% ToShift = 30 -  ToShift ;
% ToShift = -30 -  ToShift ;
% cc = circshift(cc,ToShift);

%%
% subplot(2,1,1)
% figure
line([0 0],[-(max(Rate)) max(Rate)],'linewidth',3,'color','k'); %for figure linewidth = 10
hold on
line([-(max(Rate)) max(Rate)],[0 0],'linewidth',3,'color','k');

X = smooth(x,3);
Y = smooth(y,3);
%  X1 = repmat(x,3,1);
%  Y1 = repmat(y,3,1);
% % 
% X = [];
% Y = [];
% for t= length(x):length(x)*2
%     X(t-60) = mean([X1(t-1) X1(t) X1(t+1)]);
%     Y(t-60) = mean([Y1(t-1) Y1(t) Y1(t+1)]);
% end
% 
% 
%  X = (X(1:61))';
%  Y = (Y(1:61))';

% X = (x);
% Y = (y);
hold on ;
% v = [x y];
v = [X Y];
f = (1:1:length(dirrate));
col = cc;
patch('Faces',f,'Vertices',v,'FaceVertexCData',col,'EdgeColor','interp','FaceColor','non','LineWidth',6,'LineSmoothing','off');  %for figure linewidth = 18

%for black polar plots
%  c = y;
%  patch(X,Y,c,'EdgeColor','k','LineWidth',8,'FaceColor','non');
%  line([x(n) x(1)],[y(n) y(1)],'Color',cc(n+1,:) )


hold off;
axis equal tight;
axis off;

%%

% calculate rayleigh vector length
alpha_deg = [0:binSizeDir:360-binSizeDir];
alpha_rad = deg2rad(alpha_deg)';
d=deg2rad(binSizeDir);
% [pval_angle, z_double_angle] = circ_rtest(alpha_rad, dirrate(1:60,1), d);
Vector_length = circ_r(alpha_rad, dirrate, d);
% 
%calculate PFD
PFD = circ_mean(alpha_rad, dirrate(1:60,1));
PFD_deg = rad2deg(PFD);
if PFD_deg < 0
    PFD_deg = 360+PFD_deg;

end
% 
% %peak rate
peak_HD = max(Rate);
% %mean rate
% mean_rate_HD = mean(Rate);
% 
% close gcf





% subplot(2,1,2)
% 
% Angles = [0:6:360]';
% plot(posx,posy,'Color',[.6 .6 .6],'Linewidth',2);axis square; axis off ; set(gca,'DataAspectRatio',[1 1 1]);
% 
% hold on; 
% SpikePos = [spkx spky];
% scatter(SpikePos(:,1),SpikePos(:,2));  
% 
% % SpikePos2 = posdirts;
% % for i = 1:length(SpikePos2)
% %     if SpikePos2(i) > 0
% %         SpikePos2(i) = 1;
% %     else
% %     end
% % end
% 
% % posdirts = zeros(length(posts),1);
% % posdirts(1,1) = sum(ts > 0 & ts <= posts(1,1));
% % for l = 2:size(posdirts,1)
% %     posdirts(l,1) = sum(ts > posts(l-1,1) & ts <= posts(l,1));
% % end
% 
% %dir = deg_theta(posdirts>0);
% 
% % get the direction of each spike
% dir = get_pos_spikes(ts,deg_theta,deg_theta,posts);
% 
%     
%     %scatter(OpenField.pos.xy_cm(SpikePos,1),OpenField.pos.xy_cm(SpikePos,2),'r')%dirBins = (0+anglebinwidth : anglebinwidth: 360);dirBins=dirBins'%c = linspace(1,10,360/anglebinwidth);%OpenField.pos.dir(SpikePos(1))
%     %dir = OpenField.pos.dir(SpikePos)';
%     DirbinInd = zeros(numel(dir),1);
%     nPerBin = [];
% dirBins = Angles;
%     for iBin = 1 : numel ( dirBins )
%      binInd = [];
%        if iBin ~= numel(dirBins)
%            Bin = [dirBins(iBin)  dirBins(iBin+1)];
%            if Bin(1) > Bin(2)
%                binInd = [ binInd ; find(dir > Bin(1) ); find( dir < Bin(2)) ] ;
%            else
%                binInd = [ binInd ; find(dir > Bin(1) & dir < Bin(2)   ) ] ;
%            end
%            nPerBin(iBin,1) = numel(binInd);
%            DirbinInd(binInd) = iBin;
%        else
%            Bin = [dirBins(iBin)  dirBins(1)];
%            if Bin(1) > Bin(2)
%                binInd = [ binInd ; find(dir > Bin(1) ); find( dir < Bin(2)) ] ;
%            else
%                binInd = [ binInd ; find(dir > Bin(1) & dir < Bin(2)   ) ] ;
%            end
%            nPerBin(iBin,1) = numel(binInd);
%            DirbinInd(binInd) = iBin ;
%        end
%        hold on
%        SpikesToPlot = find(DirbinInd == iBin);a = 30;
% %        scatter(posx(SpikePos(SpikesToPlot),1),posy(SpikePos(SpikesToPlot),2) ,'MarkerFaceColor', [cc(iBin,:)],'MarkerEdgeColor', [cc(iBin,:)] );
%        scatter(SpikePos(SpikesToPlot,1),SpikePos(SpikesToPlot,2) ,'MarkerFaceColor', [cc(iBin,:)],'MarkerEdgeColor', [cc(iBin,:)] );
%        hold off
%        SpikesToPlot=[];
%     end
% end

