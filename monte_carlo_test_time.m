%%
% run experiments with monte carlo  
%
%

function [mt_mean, mt_std] = monte_carlo_test_time(N, noise_sigma)

if nargin <= 0 
    N = 10; 
    noise_sigma = 1.;
end

%% generate features in the range [10, 10, 10] 
feats = createFeatures(10, 10, 10);

%% create observations for the first camera pose 
R = [1 0 0; 0 1 0; 0 0 1]; 
t = [0 0 0]'; 
cam = get_struct_core(); 
obs_i = createObservations(feats, cam, R, t); 
Pose0 = [R, t]; 

%%
gold_ea = [];
tran_ea = [];
samp_ea = []; 

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
    
    %% compute time consumption  
    tic; 
    compute_golden_error(obs, cam, R, t);
    time_gold = toc; 
    gold_ea = [gold_ea; time_gold]; 
    tic; 
    compute_transfer_error(obs, cam, R, t);
    time_tran = toc; 
    tran_ea = [tran_ea; time_tran]; 
    tic; 
    compute_sampson_error(obs, cam, R, t);
    time_sa = toc;
    samp_ea  = [time_sa; samp_ea]; 
    fprintf('monte_carlo_test_time: num: %d mean time cost: gold %f transfer: %f sampson: %f \r\n', k, time_gold, time_tran, time_sa);
end

% fprintf('final_compute_golden_error: num: %d mean: %f std: %f', size(err_array,1), mean(err_array), std(err_array));
mt_mean = [mean(gold_ea), mean(tran_ea), mean(samp_ea)]'; 
mt_std = [std(gold_ea), std(tran_ea), std(samp_ea)]';

end




