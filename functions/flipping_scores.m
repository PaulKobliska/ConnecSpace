function [Vector_length_double,flipping_score,flipping_score_range,Vector_length_quadri,Four_dir_score,Four_dir_score_minmax,Four_dir_score_range,Four_dir_score_minmax_range,correlation,minimums,maximums,corr_value,corr_bin]=flipping_scores(dirrate,binSizeDir,PFD_deg)


%% calculate double rayleigh vector length

alpha_deg_double = [0:binSizeDir*2:(360-binSizeDir)*2];
alpha_rad_double = deg2rad(alpha_deg_double)';
d_double=deg2rad(binSizeDir*2);
% [pval_angle, z_double_angle] = circ_rtest(alpha_rad, dirrate(1:60,1), d);
Vector_length_double = circ_r(alpha_rad_double, dirrate, d_double);

%% Flip score
%     subplot(2,4,4);
    circRm=dirrate;
    correlation = zeros(numel(circRm),2); %créer une 2ème colonne 
    for ibins = 0:numel(circRm)-1 
    b=circshift(circRm,ibins);%shifter de 1 bin pour pouvoir faire une corrélation en boucle
    correlation(ibins+1,1) = corr(circRm,b,'type','pearson');%calculer la corrélation entre le vecteur circRm et le vecteur shifter
    correlation(ibins+1,2) = ibins * binSizeDir; % écrire dans la deuxième colonne les angles
    end
    correlation_v2 = [correlation(:,1);1];
    correlation_v3 = [correlation(:,2);360];
    correlation = [correlation_v2 correlation_v3];

    % Calculer le flip score
    %reporter la valeur de corrélation à 180°


    correl_180 = correlation(find(correlation(:,2)==180),1);
    %reporter la valeur de corrélation à 90° (& 270 si utilise convolution)
    correl_90 = correlation(find(correlation(:,2)==90),1);
    flipping_score = correl_180-correl_90;
    









    
%     flipping_score2=sprintf('%.2f', flipping_score);

    correl_range_180 = max(correlation(29:33,1));
    correl_range_90 = max(correlation(14:18,1));
    flipping_score_range = correl_range_180-correl_range_90;
   

%% calculate quadri rayleigh vector length

alpha_deg_quadri = [0:binSizeDir*4:(360-binSizeDir)*4];
alpha_rad_quadri = deg2rad(alpha_deg_quadri)';
d_quadri=deg2rad(binSizeDir*4);
% [pval_angle, z_double_angle] = circ_rtest(alpha_rad, dirrate(1:60,1), d);
Vector_length_quadri = circ_r(alpha_rad_quadri, dirrate, d_quadri);

%% Four direction score
    correl_180 = correlation(find(correlation(:,2)==180),1);
    correl_90 = correlation(find(correlation(:,2)==90),1);
    correl_48 = correlation(find(correlation(:,2)==48),1);
    correl_138 = correlation(find(correlation(:,2)==138),1);
    Four_dir_score = mean([correl_180 correl_90])-mean([correl_48 correl_138]);   
%     Four_dir_score2=sprintf('%.2f', Four_dir_score);  
    Four_dir_score_minmax = max([correl_180 correl_90])-min([correl_48 correl_138]);   
    
    correl_range_48 = min(correlation(7:10,1));
    correl_range_138 = min(correlation(22:25,1));
    Four_dir_score_range = mean([correl_range_180 correl_range_90])-mean([correl_range_48 correl_range_138]);  
    Four_dir_score_minmax_range = max([correl_range_180 correl_range_90])-min([correl_range_48 correl_range_138]);   

%%
% % % %% calculate rayleigh sur 180° 
% % % lim_min=round(PFD_deg)-90;
% % % if lim_min < 0
% % %     lim_min = 360+lim_min;
% % % end
% % % 
% % % lim_max=round(PFD_deg)+90;
% % % if lim_max > 360
% % %     lim_max = lim_max-360;
% % % end
% % % 
% % % if lim_min>lim_max
% % %     alpha_deg_test = [lim_min-binSizeDir:-binSizeDir:lim_max];
% % %     alpha_deg_test=flip(alpha_deg_test,2);
% % % else
% % %     alpha_deg_test = [lim_min:binSizeDir:lim_max-binSizeDir];
% % % end
% % % 
% % % alpha_rad_test = deg2rad(alpha_deg_test)';
% % % 
% % % ang=[];
% % % if PFD_deg>354
% % %     PFD_deg=0;
% % % end
% % % 
% % % bin = 0:binSizeDir:360-binSizeDir;
% % % for i = 1:length(bin)-1
% % %     zbew = find(bin(i)<=PFD_deg & bin(i+1)>=PFD_deg);
% % %     if zbew ==1
% % %         ang=i;
% % %     end
% % % end
% % % freq=repmat(dirrate,3,1);
% % % fenetre=freq(ang+60-15:ang+60+14);
% % % d=deg2rad(binSizeDir);
% % % % ang=deg2rad(ang);
% % % 
% % % % [pval_angle, z_double_angle] = circ_rtest(alpha_rad, dirrate(1:60,1), d);
% % % bin_deg = [0:binSizeDir*2:360-(binSizeDir*2)];
% % % bin_rad=deg2rad(bin_deg)';
% % % Vector_length_uno = circ_r(bin_rad, fenetre, d);
% % % 
% % % [pval_uno ~] = circ_rtest(bin_rad, fenetre, d);
% % % 
% % % fenetre_oppose=freq(ang+90-15:ang+90+14);
% % % Vector_length_oppose = circ_r(bin_rad, fenetre_oppose, d);
% % % 
% % % % subplot(2,1,2)
% % % [pval_oppose ~] = circ_rtest(bin_rad, fenetre_oppose, d);

%% find number of peak in the circular autocorr

corr_value = correlation(1:33,1);
corr_bin = correlation(1:33,2);
minimums = islocalmin(corr_value,'MinSeparation',9);
maximums = islocalmax(corr_value,'MinSeparation',9);





% plot(corr_bin,corr_value,corr_bin(maximums),corr_value(maximums),'b*',corr_bin(minimums),corr_value(minimums),'r*');






















