%%
% Nov. 24, 2020, He Zhang, fuyinzh@gmail.com
% compute the sampson error based on the epipolar constraint
%

function [mean_err, std_err] = compute_sampson_error_epipolar_constrain(obs, cam, R, t)

err_array = [];
for i=1:size(obs,1)
    obs_ij = obs(i); 
    err_ij = compute_sampson_error_epipolar_instance(obs_ij, cam, R, t); 
    err_array = [err_ij; err_array]; 
    
end
    
fprintf('compute_sampson_epipolar_error: num: %d mean: %f std: %f\r\n', size(obs,1), mean(err_array), std(err_array));

mean_err = mean(err_array);
std_err = std(err_array);

end

function err = compute_sampson_error_epipolar_instance(obs, cam, R, t)
    
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
    x_true = [obs.pi.x obs.pi.y obs.pj.x obs.pj.y]'; 
    dis1 = distance(cam, x1, x_true);  
    err = dis1; 

       
end

function S = skew(x)
    
S = [0 -x(3) x(2); 
     x(3) 0 -x(1); 
     -x(2) x(1) 0]; 
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
