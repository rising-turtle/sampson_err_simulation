%%
% Nov. 25 2020, He Zhang, fuyinzh@gmail.com 
% compute the transfer distance and its error compared to ground truth  

function [mean_dis, mean_err] = compute_transfer_error(obs, cam, R, t)

err_array = [];
dis_array = []; 

for i=1:size(obs,1)
    obs_ij = obs(i); 
    [dis_ij, err_ij] = compute_transfer_error_instance(obs_ij, cam, R, t); 
    dis_array = [dis_ij; dis_array]; 
    err_array = [err_ij; err_array]; 
end
    
fprintf('compute_transfer_distance: num: %d mean_dis: %f mean_err: %f\r\n', size(obs,1), mean(dis_array), mean(err_array));

mean_dis = mean(dis_array); 
mean_err = mean(err_array);

end

function [dis, err] = compute_transfer_error_instance(obs, cam, R, t)
    
    nxi = (obs.pi_n.x-cam.cx)/cam.fx; 
    nyi = (obs.pi_n.y-cam.cy)/cam.fy; 
    % nxj = (obs.pj_n.x-cam.cx)/cam.fx; 
    % nyj = (obs.pj_n.y-cam.cy)/cam.fy; 
    d = obs.di; 
    
    pti = [nxi*d nyi*d d]'; 
    ptj = R*pti + t;
    
    ptj(1) = ptj(1)/ptj(3); 
    ptj(2) = ptj(2)/ptj(3); 
    
    %% compute transfer distance and its error 
    xj = ptj(1)*cam.fx + cam.cx; 
    yj = ptj(2)*cam.fy + cam.cy; 
    % distance |hat(X) - meas(X)|^2 
    dis = (obs.pj_n.x - xj)^2 + (obs.pj_n.y - yj)^2;
    
    % error |hat(X) - true(X)|^2 
    err = (obs.pj.x - xj)^2 + (obs.pj.y - yj)^2 + (obs.pi_n.x - obs.pi.x)^2 + (obs.pi_n.y - obs.pi.y)^2; 
    
end