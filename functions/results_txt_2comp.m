function results_txt_2comp(results_folder,inputfile,p,i,peak_rate_global,peak_rate_left,peak_rate_right,...
    peak_rate_door,corrcoeff,rotated_corrcoeff,PFD_deg,Vector_length,peak_HD,PFD_deg_left,...
    Vector_length_left,peak_HD_left,PFD_deg_right,Vector_length_right,peak_HD_right,...
    PFD_deg_door,Vector_length_door,peak_HD_door,Vector_length_double,flipping_score,...
    flipping_score_range,Vector_length_quadri,Four_dir_score,Four_dir_score_minmax,...
    Four_dir_score_range,Four_dir_score_minmax_range,number_of_corr_peaks)

if exist(strcat(results_folder,'\two_comp_results.txt'))==2 %if the file exist, 
   % open the file (a=add)
   fid=fopen(strcat(results_folder,'\two_comp_results.txt'),'a');
      % and then write results
%    for i = 1:Nunit
       fprintf(fid,'%s\t',inputfile);
       fprintf(fid,'%.0f\t',p);
       fprintf(fid,'%.0f\t',i);

       fprintf(fid,'%.4f\t',peak_rate_global);
       fprintf(fid,'%.4f\t',peak_rate_left);
       fprintf(fid,'%.4f\t',peak_rate_right);
       fprintf(fid,'%.4f\t',peak_rate_door);
       fprintf(fid,'%.4f\t',corrcoeff);
       fprintf(fid,'%.4f\t',rotated_corrcoeff);

       fprintf(fid,'%.4f\t',PFD_deg);
       fprintf(fid,'%.4f\t',Vector_length);
       fprintf(fid,'%.4f\t',peak_HD);

       fprintf(fid,'%.4f\t',PFD_deg_left);
       fprintf(fid,'%.4f\t',Vector_length_left);
       fprintf(fid,'%.4f\t',peak_HD_left);

       fprintf(fid,'%.4f\t',PFD_deg_right);
       fprintf(fid,'%.4f\t',Vector_length_right);
       fprintf(fid,'%.4f\t',peak_HD_right);

       fprintf(fid,'%.4f\t',PFD_deg_door);
       fprintf(fid,'%.4f\t',Vector_length_door);
       fprintf(fid,'%.4f\t',peak_HD_door);

       fprintf(fid,'%.4f\t',Vector_length_double);
       fprintf(fid,'%.4f\t',flipping_score);
       fprintf(fid,'%.4f\t',flipping_score_range);

       fprintf(fid,'%.4f\t',Vector_length_quadri);
       fprintf(fid,'%.4f\t',Four_dir_score);
       fprintf(fid,'%.4f\t',Four_dir_score_minmax);
       fprintf(fid,'%.4f\t',Four_dir_score_range);
       fprintf(fid,'%.4f\t',Four_dir_score_minmax_range);

       fprintf(fid,'%.0f\t',number_of_corr_peaks);
       fprintf(fid,'\n');
   
   fclose(fid);
else fid=fopen(strcat(results_folder,'\two_comp_results.txt'),'w'); %if it does not exist, open the file (w=write)
 % write header
   fprintf(fid,'%s\t','file1');
   fprintf(fid,'%s\t','tetrode');
   fprintf(fid,'%s\t','cell');
   fprintf(fid,'%s\t','peak rate map');
   fprintf(fid,'%s\t','peak rate map left');
   fprintf(fid,'%s\t','peak rate map right');
   fprintf(fid,'%s\t','peak rate map door');
   fprintf(fid,'%s\t','spatial corr');
   fprintf(fid,'%s\t','spatial corr rotated');

    fprintf(fid,'%s\t','r speed correlation');
    fprintf(fid,'%s\t','pvalue speed correlation');
    fprintf(fid,'%s\t','slope speed correlation');
    fprintf(fid,'%s\t','intercept speed correlation');

   fprintf(fid,'%s\t','PFD global');
   fprintf(fid,'%s\t','Vector length global');
   fprintf(fid,'%s\t','Peak HD glbal');
   fprintf(fid,'%s\t','PDF comp 1');
   fprintf(fid,'%s\t','Vector length comp 1');
   fprintf(fid,'%s\t','Peak HD comp 1');
   fprintf(fid,'%s\t','PDF comp 2');
   fprintf(fid,'%s\t','Vector length comp 2');
   fprintf(fid,'%s\t','Peak HD comp 2');
   fprintf(fid,'%s\t','PDF door');
   fprintf(fid,'%s\t','Vector length door');
   fprintf(fid,'%s\t','Peak HD door');
   fprintf(fid,'%s\t','Double angle vector length');
   fprintf(fid,'%s\t','Flipscore');
   fprintf(fid,'%s\t','Flipscore within range');
   fprintf(fid,'%s\t','Four angle vector length');
   fprintf(fid,'%s\t','Quadriscore');
   fprintf(fid,'%s\t','Quadriscore minmax');
   fprintf(fid,'%s\t','Quadriscore within range');
   fprintf(fid,'%s\t','Quadriscore minmax within range');
   fprintf(fid,'%s\t','Number of peaks circular autocorrelagram');

   fprintf(fid,'\n'); 
%     for i = 1:Nunit
     % add results
       fprintf(fid,'%s\t',inputfile);
       fprintf(fid,'%.0f\t',p);
       fprintf(fid,'%.0f\t',i);
       fprintf(fid,'%.4f\t',peak_rate_global);
       fprintf(fid,'%.4f\t',peak_rate_left);
       fprintf(fid,'%.4f\t',peak_rate_right);
       fprintf(fid,'%.4f\t',peak_rate_door);
       fprintf(fid,'%.4f\t',corrcoeff);
       fprintf(fid,'%.4f\t',rotated_corrcoeff);
       fprintf(fid,'%.4f\t',PFD_deg);
       fprintf(fid,'%.4f\t',Vector_length);
       fprintf(fid,'%.4f\t',peak_HD);
       fprintf(fid,'%.4f\t',PFD_deg_left);
       fprintf(fid,'%.4f\t',Vector_length_left);
       fprintf(fid,'%.4f\t',peak_HD_left);
       fprintf(fid,'%.4f\t',PFD_deg_right);
       fprintf(fid,'%.4f\t',Vector_length_right);
       fprintf(fid,'%.4f\t',peak_HD_right);
       fprintf(fid,'%.4f\t',PFD_deg_door);
       fprintf(fid,'%.4f\t',Vector_length_door);
       fprintf(fid,'%.4f\t',peak_HD_door);
       fprintf(fid,'%.4f\t',Vector_length_double);
       fprintf(fid,'%.4f\t',flipping_score);
       fprintf(fid,'%.4f\t',flipping_score_range);
       fprintf(fid,'%.4f\t',Vector_length_quadri);
       fprintf(fid,'%.4f\t',Four_dir_score);
       fprintf(fid,'%.4f\t',Four_dir_score_minmax);
       fprintf(fid,'%.4f\t',Four_dir_score_range);
       fprintf(fid,'%.4f\t',Four_dir_score_minmax_range);
       fprintf(fid,'%.0f\t',number_of_corr_peaks);

       fprintf(fid,'\n');
  fclose(fid);
end
