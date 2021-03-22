%%
% Nov.25, 2020, He Zhang, fuyinzh@gmail.com
% compute the sampson distance and its error compared to ground truth  
%

function [mean_dis, mean_err] = compute_sampson_error(obs, cam, R, t)

err_array = [];
dis_array = []; 

for i=1:size(obs,1)
    obs_ij = obs(i); 
    [dis_ij, err_ij] = compute_sampson_error_instance(obs_ij, cam, R, t); 
    dis_array = [dis_ij; dis_array]; 
    err_array = [err_ij; err_array]; 
end
    
fprintf('compute_sampson_perspective_projection: num: %d mean_dis: %f mean_err: %f\r\n', size(obs,1), mean(dis_array), mean(err_array));

mean_dis = mean(dis_array); 
mean_err = mean(err_array);

end

function [dis, err] = compute_sampson_error_instance(obs, cam, R, t)
    
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
    
    % error |hat(X) - true(X)|^2 
    x_true = [obs.pi.x obs.pi.y obs.pj.x obs.pj.y]'; 
    err = distance(cam, x1, x_true);  
    
    
  %  fprintf('sampson_error: projection error before %f after %f\n', projection_error(R,t, x, d), projection_error(R,t, x1, d)); 
    
    
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
