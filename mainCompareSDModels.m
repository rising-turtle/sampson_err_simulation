%%
% May. 15 2024, He Zhang, zhanghd@gcc.edu 
% Use sampson error to calculate pose change  
% 

function mainCompareSDModels()

% Open a file for writing
fileID = fopen('result/comparison_sd.txt', 'w');
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
   [mean_sd, std_sd] = compareSDModelsWithNoise(cam, obs_i, feats, noise_sigma);
   fprintf(fileID, '%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n', noise_sigma, mean_sd(1), std_sd(1), ...
       mean_sd(2), std_sd(2), mean_sd(3), std_sd(3), mean_sd(4), std_sd(4));
   fprintf('%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n', noise_sigma, mean_sd(1), std_sd(1), ...
       mean_sd(2), std_sd(2), mean_sd(3), std_sd(3), mean_sd(4), std_sd(4));
end


end


%% calculate mean_sd, std_sd 
function [mean_sd, std_sd] = compareSDModelsWithNoise(cam, obs_i, feats, noise_sigma)
    %% statistical data 
    sd_pro_array = []; 

    R = [1 0 0; 0 1 0; 0 0 1]; 
    t = [0 0 0]'; 
    Pose0 = [R, t]; 

    for k = 0 : 500

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
    x = estimatePoseChangeWithSD(obs, cam, false, gt); % with_epipolar = false
    x2 = estimatePoseChangeWithSD(obs, cam, true, gt); % with_epipolar = true
    % x(1:3) = x(1:3)*180.0/pi; 
    % fprintf('%f ', x); 
    
    if norm(x2) > 2
        continue;
    end

    [e_r, e_t] = calError(x, R, t); 
    [e_r2, e_t2] = calError(x2, R, t); 

    % put together 
    sd_pro_array = [e_r, e_t, e_r2, e_t2; sd_pro_array]; 

   end 
    
   % calculate mean and std 
   mean_sd = mean(sd_pro_array); 
   std_sd = std(sd_pro_array); 
    
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
function x = estimatePoseChangeWithSD(obs, cam, with_epipolar, x_true)

% x0 = [0 0 0 0.02 0.02 0.02]'; 
x0 = x_true;
% if with_epipolar
%     x0 = x_true;
% else
%     x0 = [0 0 0 0 0 0]'; 
% end
options = optimoptions('fmincon', 'Display', 'off'); % 'iter'
[x, fval] = fmincon(@(x)minSampson(x, cam, obs, with_epipolar), x0, [], [], [], [], [], [], [], options); 


% fprintf('x = %f fval = %f\n', x, fval); 
% disp(x); 

end

%% compute the sampson error function
function f = minSampson(x, cam, obs, with_epipolar)
    
f = 0; 
for i=1:size(obs,1)
    obs_ij = obs(i); 
    e = x(1:3); 
    t = x(4:6); 
    R = e2R(e);
    if with_epipolar
        dis_ij = compute_sampson_error_epipolar_instance(obs_ij, cam, R, t); 
    else
        dis_ij = compute_sampson_error_instance(obs_ij, cam, R, t); 
    end
    f = f + dis_ij; 
end


end

%% projection model based sampson_distance 
function [dis] = compute_sampson_error_instance(obs, cam, R, t)
    
    nxi = (obs.pi_n.x-cam.cx)/cam.fx; 
    nyi = (obs.pi_n.y-cam.cy)/cam.fy; 
    nxj = (obs.pj_n.x-cam.cx)/cam.fx; 
    nyj = (obs.pj_n.y-cam.cy)/cam.fy; 
    % d = obs.di; 
    d = obs.di_n;
    x = [nxi nyi nxj nyj]'; 
    e = [0 0]';
    zj = (R(3,1)*x(1)*d + R(3,2)*x(2)*d + R(3,3)*d + t(3));
    e(1) = -x(4)*zj+(R(2,1)*x(1)*d + R(2,2)*x(2)*d + R(2,3)*d + t(2));
    e(2) = x(3)*zj-(R(1,1)*x(1)*d + R(1,2)*x(2)*d + R(1,3)*d + t(1)); 
    
    de_dpi = [0 1 -nyj; -1 0 nxj]*R*[d 0; 0 d; 0 0]; 
    de_dpj = [0 -zj; zj 0]; 
    J = [de_dpi de_dpj];
    JJ_inv = (J*J');
    dx = -J'/JJ_inv*e;
    
    x1 = x + dx;
    % x2 = x - dx; 
    
    %% compute distance and error
    % distance |hat(X) - meas(X)|^2 
    x_meas = [obs.pi_n.x obs.pi_n.y obs.pj_n.x obs.pj_n.y]';
    dis = distance(cam, x1, x_meas);
    
end

function dis = distance(cam, x, xt)
    x1 = x(1) * cam.fx + cam.cx; 
    x2 = x(2) * cam.fy + cam.cy; 
    x3 = x(3) * cam.fx + cam.cx; 
    x4 = x(4) * cam.fy + cam.cy; 

    x_est = [x1 x2 x3 x4]'; 
    dis = x_est - xt; 
    dis = norm(dis)^2; 
    
end


function [dis] = compute_sampson_error_epipolar_instance(obs, cam, R, t)
    
    tji = -R' * t; 
    E = R * skew(tji); 

    nxi = (obs.pi_n.x-cam.cx)/cam.fx; 
    nyi = (obs.pi_n.y-cam.cy)/cam.fy; 
    nxj = (obs.pj_n.x-cam.cx)/cam.fx; 
    nyj = (obs.pj_n.y-cam.cy)/cam.fy; 
    
    % compute sampson distance based on epipolar constraint
    Xi = [nxi nyi 1]'; 
    Xj = [nxj nyj 1]'; 
    pXj = E'*Xj; 
    pXi = E*Xi; 
    x = [nxi nyi nxj nyj]'; 
    J = [pXj(1), pXj(2), pXi(1), pXi(2)]; 
    J2 = norm(J)^2; 
    de = Xj'*E*Xi; 
    dx = -J'*de/J2; 
    x1 = x + dx; 
    
    %% compute distance and error
    % distance |hat(X) - meas(X)|^2 
    x_meas = [obs.pi_n.x obs.pi_n.y obs.pj_n.x obs.pj_n.y]';
    dis = distance(cam, x1, x_meas); 
      
end

function S = skew(x)
    
S = [0 -x(3) x(2); 
     x(3) 0 -x(1); 
     -x(2) x(1) 0]; 
end 