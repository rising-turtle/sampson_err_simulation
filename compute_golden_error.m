%%
% Nov. 25 2020, He Zhang, fuyinzh@gmail.com 
% use golden standard way to correct the observed features 
% and compute golden standard distance and the error between ground truth 
%

function [mean_dis, mean_err] = compute_golden_error(obs, cam, R, t)

err_array = [];
dis_array = []; 

for i=1:size(obs,1)
    obs_ij = obs(i); 
    [dis_ij, err_ij] = compute_golden_error_instance(obs_ij, cam, R, t); 
    dis_array = [dis_ij; dis_array]; 
    err_array = [err_ij; err_array]; 
end
    
fprintf('compute_gold_standard: num: %d mean_dis: %f mean_err: %f\r\n', size(obs,1), mean(dis_array), mean(err_array));

mean_dis = mean(dis_array); 
mean_err = mean(err_array);

end

function f = obj_function(x, x0)
    f = x-x0;
    f = norm(f)^2;
end

function [c, ceq] = confuneq(x, d, R, t) % algebraic distance 

    % inequality constraints 
    c = []; 
    
    % equality constraints 
    
    ceq(1) = x(3)*(R(3,1)*x(1)*d + R(3,2)*x(2)*d + R(3,3)*d + t(3))-(R(1,1)*x(1)*d ...
        + R(1,2)*x(2)*d + R(1,3)*d + t(1));
    
    ceq(2) = x(4)*(R(3,1)*x(1)*d + R(3,2)*x(2)*d + R(3,3)*d + t(3))- ... 
         (R(2,1)*x(1)*d + R(2,2)*x(2)*d + R(2,3)*d + t(2));

end

function [c, ceq] = confuneq_geometric_dis(x, d, R, t) % geometric distance 

    % inequality constraints 
    c = []; 
    
    % equality constraints 
    pti = [x(1)*d, x(2)*d, d]'; 
    ptj = R * pti + t; 
    
    ceq(1) = x(3) - ptj(1)/ptj(3); 
    ceq(2) = x(4) - ptj(2)/ptj(3); 
    
    % ceq(1) = x(3)*(R(3,1)*x(1)*d + R(3,2)*x(2)*d + R(3,3)*d + t(3))-(R(1,1)*x(1)*d ...
    %    + R(1,2)*x(2)*d + R(1,3)*d + t(1));
    
    % ceq(2) = x(4)*(R(3,1)*x(1)*d + R(3,2)*x(2)*d + R(3,3)*d + t(3))- ... 
    %     (R(2,1)*x(1)*d + R(2,2)*x(2)*d + R(2,3)*d + t(2));

end

function [dis, err] = compute_golden_error_instance(obs, cam, R, t)
   
    nxi = (obs.pi_n.x-cam.cx)/cam.fx; 
    nyi = (obs.pi_n.y-cam.cy)/cam.fy; 
    nxj = (obs.pj_n.x-cam.cx)/cam.fx; 
    nyj = (obs.pj_n.y-cam.cy)/cam.fy; 
    
    x0 = [nxi nyi nxj nyj]'; 
    d = obs.di; 
  
    %% constraints 
    A = []; 
    b = []; 
    Aeq = [];
    beq = []; 
    lb = []; 
    ub = [];
    % options = optimset('Display', 'iter');
    options = optimoptions('fmincon','Display','off');  
    % [x, fval] = fmincon(@(x)obj_function(x,x0), x0, A, b, Aeq, beq, lb, ub, @(x)confuneq(x, d, R, t), options); 
    
    [x, fval] = fmincon(@(x)obj_function(x,x0), x0, A, b, Aeq, beq, lb, ub, @(x)confuneq_geometric_dis(x, d, R, t), options); 
    
    % fprintf('x = %f fval = %f\n', x, fval);
    
    %% compute distance and error 
    xi = x(1)*cam.fx + cam.cx; 
    yi = x(2)*cam.fy + cam.cy; 
    xj = x(3)*cam.fx + cam.cx; 
    yj = x(4)*cam.fy + cam.cy;
    
    % distance |hat(X) - meas(X)|^2 
    dis = (obs.pi_n.x - xi)^2 + (obs.pi_n.y - yi)^2 + (obs.pj_n.x - xj)^2 + (obs.pj_n.y - yj)^2;
    % err |hat(X) - true(X)|^2
    err = (obs.pi.x - xi)^2 + (obs.pi.y - yi)^2 + (obs.pj.x - xj)^2 + (obs.pj.y - yj)^2;
    
   % fprintf('pixel: pi (%f, %f), pi_n (%f, %f) estimate (%f, %f), \n pj (%f, %f), pj_n (%f, %f) estimate (%f, %f) err: %f\n', obs.pi.x, obs.pi.y, ...
   %   obs.pi_n.x, obs.pi_n.y, xi, yi, obs.pj.x, obs.pj.y, obs.pj_n.x, obs.pj_n.y, xj, yj, err);
    
    %% find the smallest 
%     err = 1000;
%     for i = 1:size(xi)
%         d_err = (obs.pi.x - xi)^2 + (obs.pi.y - yi)^2 + (obs.pj.x - xj)^2 + (obs.pj.y - yj)^2; 
%         if d_err < err
%             err = d_err; 
%         end
%     end
end



% 
% function err = compute_golden_error_instance(obs, cam, R, t)
%     syms pix piy pjx pjy lambda1 lambda2
%     
%     nxi = (obs.pi_n.x-cam.cx)/cam.fx; 
%     nyi = (obs.pi_n.y-cam.cy)/cam.fy; 
%     nxj = (obs.pj_n.x-cam.cx)/cam.fx; 
%     nyj = (obs.pj_n.y-cam.cy)/cam.fy; 
%     
%     F = (nxi-pix)^2 + (nyi-piy)^2 + (nxj-pjx)^2 + (nyj-pjy)^2;
%     % pj' = R*pi*d + t
%     d = obs.di; 
%     G1 = (pjx*(R(3,1)*pix*d + R(3,2)*piy*d + R(3,3)*d + t(3))-(R(1,1)*pix*d ...
%         + R(1,2)*piy*d + R(1,3)*d + t(1))) == 0;
% %     G = (pjx*(R(3,1)*pix*d + R(3,2)*piy*d + R(3,3)*d + t(3))-(R(1,1)*pix*d ...
% %         + R(1,2)*piy*d + R(1,3)*d + t(1)))^2 + ...
% %         (pjy*(R(3,1)*pix*d + R(3,2)*piy*d + R(3,3)*d + t(3))- ... 
% %         (R(2,1)*pix*d + R(2,2)*piy*d + R(2,3)*d + t(2)))^2 == 0;
%     G2 = (pjy*(R(3,1)*pix*d + R(3,2)*piy*d + R(3,3)*d + t(3))- ... 
%          (R(2,1)*pix*d + R(2,2)*piy*d + R(2,3)*d + t(2))) == 0;
%     L = F + lambda1 * lhs(G1) + lambda2 * lhs(G2);
%     dL_dpix = diff(L, pix) == 0; 
%     dL_dpiy = diff(L, piy) == 0;
%     dL_dpjx = diff(L, pjx) == 0; 
%     dL_dpjy = diff(L, pjy) == 0;
%     dL_dlambda1 = diff(L, lambda1) == 0;
%     dL_dlambda2 = diff(L, lambda2) == 0;
%     system = [dL_dpix; dL_dpiy; dL_dpjx; dL_dpjy; dL_dlambda1; dL_dlambda2]; 
%     % [xi, yi, xj, yj, la1, la2, params, conds] = solve(system, [pix piy pjx pjy lambda1 lambda2], 'ReturnConditions', true); % 'Real', true); 
%     S = solve(system, [pix piy pjx pjy lambda1 lambda2], 'MaxDegree', 3); % 'Real', true); 
% 
%     result_num = double([xi, yi, xj, yj]); 
%     
%     %% find the smallest 
%     err = 1000;
%     for i = 1:size(xi)
%         d_err = (obs.pi.x - xi)^2 + (obs.pi.y - yi)^2 + (obs.pj.x - xj)^2 + (obs.pj.y - yj)^2; 
%         if d_err < err
%             err = d_err; 
%         end
%     end
% end

% function err = compute_golden_error_instance(obs, cam, R, t)
%     syms pix piy pjx pjy lambda
%     F = (obs.pi_n.x-pix)^2 + (obs.pi_n.y-piy)^2 + (obs.pj_n.x-pjx)^2 + (obs.pj_n.y-pjy)^2;
%     % pj' = R*pi*d + t
%     d = obs.di; 
%     G = (((pjx-cam.cx)/cam.fx)*(R(3,1)*((pix-cam.cx)/cam.fx)*d + ...
%         R(3,2)*((piy-cam.cy)/cam.fy)*d + R(3,3)*d + t(3))-(R(1,1)*((pix-cam.cx)/cam.fx)*d ...
%         + R(1,2)*((piy-cam.cy)/cam.fy)*d + R(1,3)*d + t(1)))^2 + ...
%         (((pjy-cam.cy)/cam.fy)*(R(3,1)*((pix-cam.cx)/cam.fx)*d + R(3,2)*((piy-cam.cy)/cam.fy)*d + R(3,3)*d + t(3))- ... 
%         (R(2,1)*((pix-cam.cx)/cam.fx)*d + R(2,2)*((piy-cam.cy)/cam.fy)*d + R(2,3)*d + t(2)))^2 == 0;
%     L = F + lambda * lhs(G);
%     dL_dpix = diff(L, pix) == 0; 
%     dL_dpiy = diff(L, piy) == 0;
%     dL_dpjx = diff(L, pjx) == 0; 
%     dL_dpjy = diff(L, pjy) == 0;
%     dL_dlambda = diff(L, lambda) == 0;
%     system = [dL_dpix; dL_dpiy; dL_dpjx; dL_dpjy; dL_dlambda]; 
%     [xi, yi, xj, yj, la] = solve(system, [pix, piy, pjx, pjy, lambda], 'Real', true); 
%     result_num = double([xi, yi, xj, yj]); 
%     
%     %% find the smallest 
%     err = 1000;
%     for i = 1:size(xi)
%         d_err = (obs.pi.x - xi)^2 + (obs.pi.y - yi)^2 + (obs.pj.x - xj)^2 + (obs.pj.y - yj)^2; 
%         if d_err < err
%             err = d_err; 
%         end
%     end
% end