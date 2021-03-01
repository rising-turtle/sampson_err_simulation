%%
% main entrance for error metric comparison 
%
function main_error_comp()

N = 500; 
% N = 5; 
noise = 0.2:0.2:1.2; % 2.4 
% noise = 0.2:0.8:2.;
y_mean_dis = []; 
y_mean_err = [];

for i = 1:size(noise,2)
    
    sigma = noise(i);
    
    [mean_dis, mean_err] = monte_carlo_test(N, sigma); 
    y_mean_dis = [y_mean_dis; mean_dis']; 
    y_mean_err = [y_mean_err; mean_err']; 
end

%% dump to files 
% fmean = fopen('result/mean_err_gold_trans_samp.log', 'w'); 
% fstd = fopen('result/std_err_gold_trans_samp.log', 'w'); 

% fmean = fopen('result/geometric_dis_mean_err_gold_trans_samp.log', 'w'); 
% fstd = fopen('result/geometric_dis_std_err_gold_trans_samp.log', 'w'); 

% fmean = fopen('result/sampson_dis_epipolar_vs_perspective.log', 'w'); 

%% update gold standard 
% fmean_dis = fopen('result/mean_residual_error_comparison.log', 'w'); 
% fmean_err = fopen('result/mean_estimation_error_comparison.log', 'w'); 

fmean_dis = fopen('result/mean_residual_error_comparison_new.log', 'w'); 
fmean_err = fopen('result/mean_estimation_error_comparison_new.log', 'w'); 

fprintf(fmean_dis, '%3.3f  %7.7f  %7.7f  %7.7f %7.7f \r\n', [noise', y_mean_dis(:,1), y_mean_dis(:,2), y_mean_dis(:,3), y_mean_dis(:,4)]');
fprintf(fmean_err, '%3.3f  %7.7f  %7.7f  %7.7f %7.7f \r\n', [noise', y_mean_err(:,1), y_mean_err(:,2), y_mean_err(:,3), y_mean_err(:,4)]');

% fprintf(fmean, '%3.3f  %7.7f  %7.7f \r\n', [noise', y_mean(:,1), y_mean(:,2)]');
% fprintf(fstd, '%3.3f  %7.7f  %7.7f  %7.7f %7.7f \r\n', [noise', y_std(:,1), y_std(:,2), y_std(:,3), y_std(:,4)]');

fclose(fmean_dis); 
fclose(fmean_err); 

end