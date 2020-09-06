%%
% main entrance for error metric comparison 
%
function main_error_comp()

N = 500; 
noise = 0.2:0.2:5;
y_mean = []; 
y_std = [];

for i = 1:size(noise,2)
    
    sigma = noise(i);
    
    [e_mean, e_std] = monte_carlo_test(N, sigma); 
    y_mean = [y_mean; e_mean']; 
    y_std = [y_std; e_std']; 
end

%% dump to files 
% fmean = fopen('result/mean_err_gold_trans_samp.log', 'w'); 
% fstd = fopen('result/std_err_gold_trans_samp.log', 'w'); 

fmean = fopen('result/geometric_dis_mean_err_gold_trans_samp.log', 'w'); 
fstd = fopen('result/geometric_dis_std_err_gold_trans_samp.log', 'w'); 

fprintf(fmean, '%3.3f  %7.7f  %7.7f  %7.7f %7.7f \r\n', [noise', y_mean(:,1), y_mean(:,2), y_mean(:,3), y_mean(:,4)]');
fprintf(fstd, '%3.3f  %7.7f  %7.7f  %7.7f %7.7f \r\n', [noise', y_std(:,1), y_std(:,2), y_std(:,3), y_std(:,4)]');

fclose(fmean); 
fclose(fstd); 

end