%%
% Nov. 25 2020, He Zhang, fuyinzh@gmail.com 
% run experiments with monte carlo  
%
%

function [mt_mean_dis, mt_mean_err] = monte_carlo_test(N, noise_sigma)

if nargin <= 0 
    N = 10; 
    noise_sigma = 1.;
end

%% generate features in the range [10, 10, 10] 
feats = createFeatures(10, 10, 10);

%% create observations for the first camera pose 
R = [1 0 0; 0 1 0; 0 0 1]; 
t = [0 0 0]'; 
Pose0 = [R, t]; 
cam = get_struct_core(); 
obs_i = createObservations(feats, cam, R, t); 

%%
re_dis_array = []; re_err_array = [];
td_dis_array = []; td_err_array = []; 
sd_dis_array = []; sd_err_array = [];  
sd_epi_dis_array = []; sd_epi_err_array = []; 


%% random 
for k =0:N
    % randomly generate next camera pose within rpy [-30, 30] degree, xyz
    rr = randn(6,1); 
    s = max(rr); 
    if s < -min(rr)
        s = -min(rr);
    end
    rr(1:3) = rr(1:3)/s*30; 
    rr(4:6) = rr(4:6)/s*2; 
    euler_angle = rr(1:3)*pi/180;
    R = e2R(euler_angle);
    t = rr(4:6); 
    Pose1 = [R, t];
    %% find observations 
    obs_j = createObservations(feats, cam, R, t);
    
    %% overlap observations
    obs = overlap_obs(obs_i, obs_j); 
    if size(obs,1) < 10
        fprintf('too small observations for k = %f with %d obs', k+1, size(obs,1));
        continue; 
    end
    
    %% add noise 
    obs = add_noise(obs, noise_sigma); 
    
    %% add normalized feature measurement 
    obs = normalize_feature_measure(obs, cam); 
    
    %% add triangulated depth 
    obs = triangulate_depth(Pose0, Pose1, obs); 
    
    %% compute error 
    [re_dis, re_err] = compute_golden_error(obs, cam, R, t);
    [td_dis, td_err] = compute_transfer_error(obs, cam, R, t);
    [sd_dis, sd_err] = compute_sampson_error(obs, cam, R, t);
    % [sg_me, ] = compute_sampson_error_geometric_dis(obs, cam, R, t); 
    [sd_epi_dis, sd_epi_err] = compute_sampson_error_epipolar_constrain(obs, cam, R, t); 
    
    re_dis_array = [re_dis; re_dis_array]; re_err_array = [re_err; re_err_array]; 
    td_dis_array = [td_dis; td_dis_array]; td_err_array = [td_err; td_err_array];  
    sd_dis_array = [sd_dis; sd_dis_array]; sd_err_array = [sd_err; sd_err_array]; 
    sd_epi_dis_array = [sd_epi_dis; sd_epi_dis_array]; sd_epi_err_array = [sd_epi_err; sd_epi_err_array]; 

end

% fprintf('final_compute_golden_error: num: %d mean: %f std: %f', size(err_array,1), mean(err_array), std(err_array));
mt_mean_dis = [mean(re_dis_array), mean(td_dis_array), mean(sd_dis_array), mean(sd_epi_dis_array)]'; 
mt_mean_err = [mean(re_err_array), mean(td_err_array), mean(sd_err_array), mean(sd_epi_err_array)]';

% mt_mean_dis = [mean(samp_ea), mean(samp_epi_ea)]'; 
% mt_std = [std(samp_ea), std(samp_epi_ea)]'; 

end




