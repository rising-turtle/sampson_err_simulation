%%
% May. 15 2024, He Zhang, zhanghd@gcc.edu 
% Use transfer error to calculate pose change  
% 

function mainCompareTD()

% Open a file for writing
fileID = fopen('result/comparison_td.txt', 'w');
% Check if the file was opened successfully
if fileID == -1
    error('Failed to open the file.');
    return; 
end

%% image noise 
noise_sigma = 1.; % 1pixel 

%% generate features in the range [10, 10, 10] 
% feats = createFeatures(10, 10, 10);
feats = createFeatures(5, 5, 5);

%% create observations for the first camera pose 
cam = get_struct_core(); 
R = [1 0 0; 0 1 0; 0 0 1]; 
t = [0 0 0]'; 
obs_i = createObservations(feats, cam, R, t); 

noise = 0.2:0.2:3.0; % 2.4 

for i = 1:size(noise,2)
   noise_sigma = noise(i);
   [mean_td, std_td] = compareTDModelsWithNoise(cam, obs_i, feats, noise_sigma);
   fprintf(fileID, '%f\t%f\t%f\t%f\t%f\n', noise_sigma, mean_td(1), std_td(1), ...
       mean_td(2), std_td(2));
   fprintf('%f\t%f\t%f\t%f\t%f\n', noise_sigma, mean_td(1), std_td(1), ...
       mean_td(2), std_td(2));
end


end


%% calculate mean_sd, std_sd 
function [mean_td, std_td] = compareTDModelsWithNoise(cam, obs_i, feats, noise_sigma)
    %% statistical data 
    td_pro_array = []; 

    R = [1 0 0; 0 1 0; 0 0 1]; 
    t = [0 0 0]'; 
    Pose0 = [R, t]; 

    for k = 0 : 1000

    %% fake a pose change 
    rr = [5, 7, -4, 0.1, 0.2, -0.3]; 
    rr(1:3) = 5.*randn(1,3); 
    rr(4:6) = 0.2*randn(1,3);

    euler_angle = rr(1:3)*pi/180;
    R = e2R(euler_angle);
    t = rr(4:6)'; 
    Pose1 = [R, t];
    gt = [euler_angle, rr(4:6)]';
    
    %% find observations 
    obs_j = createObservations(feats, cam, R, t);
    
    %% overlap observations
    obs = overlap_obs(obs_i, obs_j); 
    if size(obs,1) < 10
        fprintf('too small observations for k = %f with %d obs', k+1, size(obs,1));
        return; 
    end
    %% add noise 
   
    obs = add_noise(obs, noise_sigma); 
    
    %% add normalized feature measurement 
    obs = normalize_feature_measure(obs, cam); 
    
    %% add triangulated depth 
    obs = triangulate_depth(Pose0, Pose1, obs); 
    
    %% compute estimate x 
    x = estimatePoseChangeWithTD(obs, cam, gt); % with_epipolar = false

    [e_r, e_t] = calError(x, R, t); 

    % put together 
    td_pro_array = [e_r, e_t; td_pro_array]; 

   end 
    
   % calculate mean and std 
   mean_td = mean(td_pro_array); 
   std_td = std(td_pro_array); 
    
end

%% calculate error 
function [e_r, e_t] = calError(x_est, R, t)
        %% calculate error 
    % http://www.boris-belousov.net/2016/12/01/quat-dist/
    R_est = e2R(x_est(1:3)); 
    R_diff = R*R_est'; 
    % Calculate the trace of the relative rotation matrix
    trace_value = trace(R_diff);

    % Compute the angle in radians between the two rotations
    e_r = acos((trace_value - 1) / 2);
    
    % err_translation 
    t_diff = x_est(4:6) - t; 
    e_t = sqrt(dot(t_diff,t_diff'));
end

%% estimate pose change given the matched features 
function x = estimatePoseChangeWithTD(obs, cam, x_true)

% x0 = [0 0 0 0.02 0.02 0.02]'; 
x0 = x_true;
options = optimoptions('fmincon', 'Display', 'off'); % 'iter'
[x, fval] = fmincon(@(x)minTD(x, cam, obs), x0, [], [], [], [], [], [], [], options); 


% fprintf('x = %f fval = %f\n', x, fval); 
% disp(x); 

end

%% compute the sampson error function
function f = minTD(x, cam, obs)
    
f = 0; 
for i=1:size(obs,1)
    obs_ij = obs(i); 
    e = x(1:3); 
    t = x(4:6); 
    R = e2R(e);
    dis_ij = compute_transfer_error_instance(obs_ij, cam, R, t); 
    f = f + dis_ij; 
end


end

%% projection model based sampson_distance 
function [dis] = compute_transfer_error_instance(obs, cam, R, t)
    
    nxi = (obs.pi_n.x-cam.cx)/cam.fx; 
    nyi = (obs.pi_n.y-cam.cy)/cam.fy; 
    % nxj = (obs.pj_n.x-cam.cx)/cam.fx; 
    % nyj = (obs.pj_n.y-cam.cy)/cam.fy; 
    % d = obs.di; 
    d = obs.di_n;
    
    pti = [nxi*d nyi*d d]'; 
    ptj = R*pti + t;
    
    ptj(1) = ptj(1)/ptj(3); 
    ptj(2) = ptj(2)/ptj(3); 
    
    %% compute transfer distance and its error 
    xj = ptj(1)*cam.fx + cam.cx; 
    yj = ptj(2)*cam.fy + cam.cy; 
    % distance |hat(X) - meas(X)|^2 
    dis = (obs.pj_n.x - xj)^2 + (obs.pj_n.y - yj)^2;
       
end

