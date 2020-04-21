%%
% given camera pose, and feature locations, generate observations
% compute feature observations given the camera pose [R,t] and 
%
function obs = createObservations(feats, cam, R_cam, t_cam)

% get feature points
gpc = get_feature_pc(feats); 

N = size(gpc,2);
% camera model 
% cam = get_struct_core(); 

obs = [];

% Rc2o = R_cam';
% t = -Rc2o * t_cam;
% lpc = transform_pc(Rc2o, t, gpc);

lpc = transform_pc(R_cam, t_cam, gpc);

for i=1:N
    %% check out how many features are available 
    fpt = lpc(:, i);
    [good, px, py] = in_cam_view(cam, fpt);
    if good == 1
        obs_ij = generate_obs(gpc(:,i), lpc(:,i), feats(i).id, px, py);
        obs = [obs; obs_ij];
    end
end

end

%% generate observation for this point
function obs = generate_obs(gpt, lpt, feat_id, px, py)
    obs.gpt = gpt;
    obs.lpt = lpt;
    % obs.pose_id = pose_i;
    obs.feat_id = feat_id;
    obs.obs_x = px;
    obs.obs_y = py; 
end

%% get point cloud from features
function gpc = get_feature_pc(feats)
    n = length(feats);
    gpc = zeros(3, n);
    for i = 1:n
        ft = feats(i); 
        gpc(1, i) = ft.x;
        gpc(2, i) = ft.y; 
        gpc(3, i) = ft.z;
    end
end

%% transform point cloud from global coordinate into camera's coordinate 
function lpc = transform_pc(R, t, gpc)
    [m, n] = size(gpc); 
    translation = repmat(t, 1, n); 
    lpc = R * gpc + translation; 
end
