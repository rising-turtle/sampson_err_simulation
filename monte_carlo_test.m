%%
% run experiments with monte carlo  
%
%

function [mt_mean, mt_std] = monte_carlo_test(N, noise_sigma)

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
    
    %% compute error 
    [g_me, ] = compute_golden_error(obs, cam, R, t);
    [t_me, ] = compute_transfer_error(obs, cam, R, t);
    [s_me, ] = compute_sampson_error(obs, cam, R, t);
    gold_ea = [g_me; gold_ea]; 
    tran_ea = [t_me; tran_ea]; 
    samp_ea = [s_me; samp_ea]; 
    
end

% fprintf('final_compute_golden_error: num: %d mean: %f std: %f', size(err_array,1), mean(err_array), std(err_array));
mt_mean = [mean(gold_ea), mean(tran_ea), mean(samp_ea)]'; 
mt_std = [std(gold_ea), std(tran_ea), std(samp_ea)]';

end

function obs = overlap_obs(obs_i, obs_j)
Ni = size(obs_i, 1); 
Nj = size(obs_j, 1);
obs = []; 
for i =1:Ni
    id_i = obs_i(i).feat_id;
    good = 0; 
    jj = -1;
    for j = 1:Nj
        id_j = obs_j(j).feat_id; 
        if id_i == id_j
            jj = j; 
            good = 1; 
            break; 
        end
    end
    
    if good == 1 % find one 
        obs_ij.pi.x = obs_i(i).obs_x; 
        obs_ij.pi.y = obs_i(i).obs_y; 
        obs_ij.di = obs_i(i).lpt(3); %% depth 
        obs_ij.pj.x = obs_j(jj).obs_x; 
        obs_ij.pj.y = obs_j(jj).obs_y; 
        obs_ij.pi_n = obs_ij.pi; 
        obs_ij.pj_n = obs_ij.pj;
        obs = [obs_ij; obs];
    end
    
end

end

function obs_n = add_noise(obs, sigma)

    obs_n = obs; 
    for i = 1:size(obs,1)
        obs_n(i).pi_n.x = obs_n(i).pi_n.x + randn(1,1)*sigma; 
        obs_n(i).pi_n.y = obs_n(i).pi_n.y + randn(1,1)*sigma; 
        obs_n(i).pj_n.x = obs_n(i).pj_n.x + randn(1,1)*sigma; 
        obs_n(i).pj_n.y = obs_n(i).pj_n.y + randn(1,1)*sigma; 
    end

end



