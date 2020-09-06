%%
% Sep. 5 2020, He Zhang, fuyinzh@gmail.com
% compute the sampson error with the geometric distance 
% and compute the error 
%

function [mean_err, std_err] = compute_sampson_error_geometric_dis(obs, cam, R, t)

err_array = [];
for i=1:size(obs,1)
    obs_ij = obs(i); 
    err_ij = compute_sampson_error_instance(obs_ij, cam, R, t); 
    err_array = [err_ij; err_array]; 
    
end
    
fprintf('compute_sampson_error_geo: num: %d mean: %f std: %f\r\n', size(obs,1), mean(err_array), std(err_array));

mean_err = mean(err_array);
std_err = std(err_array);

end

function err = compute_sampson_error_instance(obs, cam, R, t)
    
    nxi = (obs.pi_n.x-cam.cx)/cam.fx; 
    nyi = (obs.pi_n.y-cam.cy)/cam.fy; 
    nxj = (obs.pj_n.x-cam.cx)/cam.fx; 
    nyj = (obs.pj_n.y-cam.cy)/cam.fy; 
    d = obs.di; 
    x = [nxi nyi nxj nyj]'; 
    e = [0 0]';
    % zj = (R(3,1)*x(1)*d + R(3,2)*x(2)*d + R(3,3)*d + t(3));
    pti = [nxi*d, nyi*d, d]'; 
    ptj = R*pti + t;
    
    % e(1) = -x(4)*zj+(R(2,1)*x(1)*d + R(2,2)*x(2)*d + R(2,3)*d + t(2));
    % e(2) = x(3)*zj-(R(1,1)*x(1)*d + R(1,2)*x(2)*d + R(1,3)*d + t(1)); 

    e(1) = ptj(1)/ptj(3) - x(3); 
    e(2) = ptj(2)/ptj(3) - x(4); 
    
    % de_dpi = [0 1 -nyj; -1 0 nxj]*R*[d 0; 0 d; 0 0]; 
    % de_dpj = [0 -zj; zj 0]; 
    
    de_dpi = [1/ptj(3) 0 -ptj(1)/ptj(3)^2; 0 1/ptj(3) -ptj(2)/ptj(3)^2]*R*[d 0; 0 d; 0 0]; 
    de_dpj = [-1 0; 0 -1]; 
    
    J = [de_dpi de_dpj];
    JJ_inv = (J*J');
    dx = -J'/JJ_inv*e;
    
    x1 = x + dx;
    % x2 = x - dx; 
    
    %% find out which one is more close 
    x_true = [obs.pi.x obs.pi.y obs.pj.x obs.pj.y]'; 
    
    dis1 = distance(cam, x1, x_true); 
    % dis2 = distance(cam, x2, x_true); 
    
    err = dis1; % min(dis1, dis2); 
    
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
