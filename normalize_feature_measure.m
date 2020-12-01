%%
% Nov. 30 2020, He Zhang, fuyinzh@gmail.com 
% 
% normalize visual feature measurement 
%

function obs_n = normalize_feature_measure(obs, cam)
    obs_n = obs; 
    for i = 1:size(obs,1)
        obs_n(i).ni.x = norm_x(obs_n(i).pi.x, cam);
        obs_n(i).ni_n.x = norm_x(obs_n(i).pi_n.x, cam);
        
        obs_n(i).ni.y = norm_y(obs(i).pi.y, cam);
        obs_n(i).ni_n.y = norm_y(obs(i).pi_n.y, cam);
            
        obs_n(i).nj.x = norm_x(obs_n(i).pj.x, cam);
        obs_n(i).nj_n.x = norm_x(obs_n(i).pj_n.x, cam);
        
        obs_n(i).nj.y = norm_y(obs(i).pj.y, cam);
        obs_n(i).nj_n.y = norm_y(obs(i).pj_n.y, cam);
    end
end

function nx = norm_x(x, cam)
      nx = (x-cam.cx)/cam.fx; 
end

function ny = norm_y(y, cam)
      ny = (y-cam.cy)/cam.fy; 
end