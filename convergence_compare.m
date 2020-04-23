%%
% compare the distance between corrected point and true point
% 
%

function convergence_compare(N, sig_array)
    if nargin <= 0 
        N = 10; 
        sig_array = [1:0.2:4]';
    end
    
    mean_sampson_e = []; 
    mean_transfer_e = []; 
    std_sampson_e = []; 
    std_transfer_e = []; 
    
    for i = 1:size(sig_array,1)
        
        [dis_sampson, dis_transfer] = monte_carlo(N, sig_array(i)); 
        
        mean_sampson_e = [mean_sampson_e; mean(dis_sampson)];
        mean_transfer_e = [mean_transfer_e; mean(dis_transfer)]; 
        std_sampson_e = [std_sampson_e; std(dis_sampson)]; 
        std_transfer_e = [std_transfer_e; std(dis_transfer)]; 
        
        fprintf('mean_sampson_e = %f mean_transfer_e = %f noise_sigma = %f \n', ...
            mean(dis_sampson), mean(dis_transfer), sig_array(i));
%         fprintf('std_sampson_e = %f std_transfer_e = %f noise_sigma = %f \n', ...
%             std(dis_sampson), std(dis_transfer), sig_array(i));
    end
    
    %% TODO save the result into file

end

function [dis_sampson, dis_transfer] = monte_carlo(N, noise_sigma)

%% generate features in a range
feats = createFeatures(8, 8, 7);

%% create observations for the first camera pose 
R = [1 0 0; 0 1 0; 0 0 1]; 
t = [0 0 0]'; 
cam = get_struct_core(); 
obs_i = createObservations(feats, cam, R, t); 

%% number of trials that pass chi2 test 
dis_sampson = []; 
dis_transfer = []; 

%% random 
for k = 1:N
    
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
        %fprintf('too small observations for k = %f with %d obs\n', k+1, size(obs,1));
        continue; 
    end
    
    %% add noise
    obs = add_noise(obs, noise_sigma); 
    
    for i=1:10
       sampson_e = sampson_corrected_dis(obs(i), cam, R, t); 
       transfer_e = transfer_error_dis(obs(i), cam, R, t);
       dis_sampson = [dis_sampson; sampson_e]; 
       dis_transfer = [dis_transfer; transfer_e]; 
    end
end

end

function dis = transfer_error_dis(obs, cam, R, t)

    nxi = (obs.pi_n.x-cam.cx)/cam.fx; 
    nyi = (obs.pi_n.y-cam.cy)/cam.fy; 
    % nxj = (obs.pj_n.x-cam.cx)/cam.fx; 
    % nyj = (obs.pj_n.y-cam.cy)/cam.fy; 
    d = obs.di; 
    
    pti = [nxi*d nyi*d d]'; 
    ptj = R*pti + t;
    
    ptj(1) = ptj(1)/ptj(3); 
    ptj(2) = ptj(2)/ptj(3); 
    
    % found point 
    x = [nxi nyi ptj(1) ptj(2)]'; 
    x_true = [obs.pi.x obs.pi.y obs.pj.x obs.pj.y]'; 
    dis = distance(cam, x, x_true); 
    
end

function dis = sampson_corrected_dis(obs, cam, R, t)

    %% compute error
    nxi = (obs.pi_n.x-cam.cx)/cam.fx; 
    nyi = (obs.pi_n.y-cam.cy)/cam.fy; 
    nxj = (obs.pj_n.x-cam.cx)/cam.fx; 
    nyj = (obs.pj_n.y-cam.cy)/cam.fy; 
    d = obs.di; 
    x = [nxi nyi nxj nyj]'; 
    e = [0 0]';
    zj = (R(3,1)*x(1)*d + R(3,2)*x(2)*d + R(3,3)*d + t(3));
    e(1) = -x(4)*zj+(R(2,1)*x(1)*d + R(2,2)*x(2)*d + R(2,3)*d + t(2));
    e(2) = x(3)*zj-(R(1,1)*x(1)*d + R(1,2)*x(2)*d + R(1,3)*d + t(1)); 

    %% compute covariance matrix 
    de_dpi = [0 1 -nyj; -1 0 nxj]*R*[d 0; 0 d; 0 0]; 
    de_dpj = [0 -zj; zj 0]; 
    J = [de_dpi de_dpj];
    
    JJ_inv = (J*J');
    dx = -J'/JJ_inv*e;
    
    x1 = x + dx;
    
    %% find out which one is more close 
    x_true = [obs.pi.x obs.pi.y obs.pj.x obs.pj.y]'; 
    
    dis = distance(cam, x1, x_true); 
end

function dis = distance(cam, x, xt)
    x1 = x(1) * cam.fx + cam.cx; 
    x2 = x(2) * cam.fy + cam.cy; 
    x3 = x(3) * cam.fx + cam.cx; 
    x4 = x(4) * cam.fy + cam.cy; 

    x_est = [x1 x2 x3 x4]'; 
    dis = x_est - xt; 
    dis = norm(dis); 
    
end
