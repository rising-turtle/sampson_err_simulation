%%
% Nov. 30 2020, He Zhang, fuyinzh@gmail.com 
% 
% triangulate depth 
%

function obs_n = triangulate_depth(Pose0, Pose1, obs)
    obs_n = obs; 
    for i = 1:size(obs,1)
        
        pt0 = [obs(i).ni_n.x, obs(i).ni_n.y]'; 
        pt1 = [obs(i).nj_n.x, obs(i).nj_n.y]'; 
        pt3d = triangulate_point(Pose0, Pose1, pt0, pt1); 
        
        obs_n(i).di_n = pt3d(3); 
    
    end
end