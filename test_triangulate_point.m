%%
% Nov. 30 2020, He Zhang, fuyinzh@gmail.com 
% 
% test triangulate point 
%

%% generate features in the range [10, 10, 10] 
feats = createFeatures(10, 10, 10);

%% create observations for the first camera pose 
R = [1 0 0; 0 1 0; 0 0 1]; 
t = [0 0 0]'; 
Pose0 = [R, t];
cam = get_struct_core(); 
obs_i = createObservations(feats, cam, R, t); 

%% create observations for the second camera pose 
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
    fprintf('too small observations with %d obs', size(obs,1));
end

%% add noise
noise_sigma = 1.;
obs = add_noise(obs, noise_sigma);     

for i=1:size(obs,1)
    nxi = (obs(i).pi_n.x-cam.cx)/cam.fx; 
    nyi = (obs(i).pi_n.y-cam.cy)/cam.fy; 
    
    nxj = (obs(i).pj_n.x-cam.cx)/cam.fx;
    nyj = (obs(i).pj_n.y-cam.cy)/cam.fy; 
    
    pt0 = [nxi, nyi]'; 
    pt1 = [nxj, nyj]'; 
    pt3d = triangulate_point(Pose0, Pose1, pt0, pt1); 
    fprintf('feature pair i=%d true depth = %f, triangulate depth = %f \n',  ...
                i, obs(i).di, pt3d(3));
    if i > 20
        break;
    end
end









