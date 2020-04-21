%%
% compute the transfer error  
% and compute the error 
%

function [mean_err, std_err] = compute_transfer_error(obs, cam, R, t)

err_array = [];
for i=1:size(obs,1)
    obs_ij = obs(i); 
    err_ij = compute_transfer_error_instance(obs_ij, cam, R, t); 
    err_array = [err_ij; err_array]; 
    
end
    
fprintf('compute_transfer_error: num: %d mean: %f std: %f', size(obs,1), mean(err_array), std(err_array));

mean_err = mean(err_array);
std_err = std(err_array);

end

function err = compute_transfer_error_instance(obs, cam, R, t)
    
    nxi = (obs.pi_n.x-cam.cx)/cam.fx; 
    nyi = (obs.pi_n.y-cam.cy)/cam.fy; 
    % nxj = (obs.pj_n.x-cam.cx)/cam.fx; 
    % nyj = (obs.pj_n.y-cam.cy)/cam.fy; 
    d = obs.di; 
    
    pti = [nxi*d nyi*d d]'; 
    ptj = R*pti + t;
    
    ptj(1) = ptj(1)/ptj(3); 
    ptj(2) = ptj(2)/ptj(3); 
    
     %% compute the error 
    xj = ptj(1)*cam.fx + cam.cx; 
    yj = ptj(2)*cam.fy + cam.cy; 
    err = (obs.pj_n.x - xj)^2 + (obs.pj_n.y - yj)^2;
    
end