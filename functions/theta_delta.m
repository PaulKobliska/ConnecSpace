function [peak_frequency_theta,theta_power,delta_power,sig_threshold,theta_modulation] = theta_delta(TimeWindow,Tbin,ts)

time_window = TimeWindow; %longueur de l'autocorr
% Tbin = 10; %bin autocorr
Nbins1 = ceil(ts(end));
hist_spk = hist(ts, 0.5:1:Nbins1-0.5);
time_corr = xcorr(hist_spk, time_window);
% xlim([1 500]);
%%
%pour changer le xcorr default bin value
time_corr(time_window +1) = []; %enlever la valeur de la correlation max (au lag 0)
maxbin =ceil(length(time_corr)/Tbin);
subs = repmat(1:maxbin,Tbin,1);
subs = reshape(subs,Tbin*maxbin,1);
subs=subs(1:length(time_corr)); %les 2 vecteurs doivent faire la mm taille pour accumarray
time_corr = accumarray(subs,time_corr)'; %somme les valeurs de time_corr sur Tbin
newlags=floor(time_window/Tbin);

% subplot(3,4,10);
% y = -time_window:Tbin:time_window;
% y(y==0)=[];
% % figure;
% bar(y,time_corr);
% set(get(gca,'Children'),'FaceColor',[0.6 0.6 0.6]);
% set(get(gca,'Children'),'EdgeColor',[0 0 0]);
% set(gca,'TickDir','out','Box','off');
% set(get(gca,'XLabel'),'Fontsize',12,'String','Interval (ms)');
% set(get(gca,'YLabel'),'Fontsize',12,'String','N spikes');
% set(gca,'ylim',[-1,max(time_corr)*1.2]);
% %xlim([0 (TimeWindow + (TimeWindow*0.1)) ])
% xlim([-(time_window+(time_window*0.1)) (time_window + (time_window*0.1)) ])
% ylims = get(gca,'ylim');
% ylim([0 ylims(2)]);
% line([0 0],[0 ylims(2)],'Color',[0 0 0]);
% axis square
% title('Autocorrelogram');

%%
%fourier transform
Fs = (time_window / Tbin); % Sampling frequency                    
T = 1/Fs;             % Sampling period       
% L = (time_window / Tbin)*2 ;             % Length of signal
L = length(time_corr);
t = (0:L-1)*T;        % Time vector

%Compute the two-sided spectrum P2. Then compute the single-sided spectrum P1 based on P2 and the even-valued signal length L.
Y = fft(time_corr);
P2 = abs(Y/L);
P1 = P2(1:L/2+1); %single sided spectrum
P1(2:end-1) = 2*P1(2:end-1); %

%Define the frequency domain f and plot the single-sided amplitude spectrum P1.
f = Fs*(0:(L/2))/L;

 min_theta = min(find(f>= 4));
 max_theta = max(find(f<=12));
 min_delta = min(find(f>=2));
 max_delta = max(find(f<=4));

[maxt, indt] = max(P1(1,min_theta:max_theta));
peak_frequency_theta = f(min_theta+(indt-1));
theta_power = maxt;
[maxd, indd] = max(P1(1,min_delta:max_delta));
%for delta band, use local maxima function
% [maxd,indd] = findpeaks(P1(1,min_delta:max_delta));
delta_power = maxd;

if isempty(maxd)
    maxd = 0;
    indd = 0;
    peak_frequency_delta = 0;
else
    peak_frequency_delta = f(min_delta+(indd-1));
end

mean_without_theta = mean([P1(2:indt+min_theta-2) P1(indt+min_theta:end)]);

if maxd == 0
    mean_without_delta = 0;
else
    mean_without_delta = mean([P1(2:indd+min_delta-2) P1(indd+min_delta:end)]);
end

if maxd == 0
   mean_without_theta_delta = mean_without_theta;
else
   mean_without_theta_delta = mean([P1(3:indd+min_delta-2) P1(indd+min_delta:indt+min_theta-2) P1(indt+min_theta:end)]);
end

sig_threshold = 5*mean_without_theta_delta;

if theta_power>sig_threshold;
    theta_modulation = 1;
else
    theta_modulation = 0;
end    

% subplot(3,4,11);
plot(f,P1,'-k','linewidth',2,'markerfacecolor',[.7 .7 .7])
line([0 50],[sig_threshold sig_threshold],'color','r','linewidth',2);
xlim([1 50]);
set(gca,'ylim',[-1,max(P1(2:end))*1.2]);
box off
set(gca,'TickDir','out','Box','off');
set(get(gca,'XLabel'),'Fontsize',12,'String','frequency (hz)');
set(get(gca,'YLabel'),'Fontsize',12,'String','power');
axis square
title('FFT');



% %to write modulation on the figure
% theta_peak2 = sprintf('%.0f',peak_frequency_theta);
% theta_power2 = sprintf('%.2f',maxt);
% delta_peak2 = sprintf('%.0f',peak_frequency_delta);
% delta_power2 = sprintf('%.2f',maxd);
% mean_band = sprintf('%.2f',mean_without_theta_delta);
% subplot(1,3,3);
%     text(.1,.9,strcat('theta peak=',theta_peak2,'hz'));
%     text(.1,.7,strcat('theta power=',theta_power2));
%     text(.1,.5,strcat('delta peak=',delta_peak2,'hz'));    
%     text(.1,.3,strcat('delta power=',delta_power2));
%     text(.1,.1,strcat('mean band=',mean_band));
% axis off
axis square
xlim([f(min_theta) f(max_theta)]);

end