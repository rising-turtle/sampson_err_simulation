%%
%   April, 30, 2020, He Zhang, hzhang8@vcu.edu 
% compute jacobian of the matrix JJt 
%       where J =[de_dpti; de_dptj] 
%              de_dpti = [1 0 -ptj.y; -1 0 ptj.x] * Rj'* Ri * [d 0; 0 d; 0 0]
%              de_dptj = [0 -zj; zj 0]
%              

function [dJ_dpi, dJ_dpj] = jacobian_JJt(pose_i, pose_j, d, pti, ptj)

if nargin <= 0
    pose_i = [1 2 3 sqrt(1-0.1^2-0.3^2-0.5^2) 0.1 -0.3 0.5]'; 
    pose_j = [0.1 0.2 0.3 sqrt(1-0.4^2-0.1^2-0.2^2) 0.4 0.1 -0.2]';
    d = 3;
    pti = [0.1, 0.2]'; 
    ptj = [-0.3, 1.]';
end
    J = J_pose(pose_i, pose_j, d, pti, ptj);
    Jt = J'; 
    [dJ_dpose_i, dJ_dpose_j] = jacobian_J(pose_i, pose_j, d, pti, ptj);
    [dJt_dpose_i, dJt_dpose_j] = jacobian_J_transpose(pose_i, pose_j, d, pti, ptj);
    
    %% numeric result 
    dJJt_d_pose_i_num = numeric_jacobian_matrix_pose(@(pose_i) JJt(pose_i, pose_j, d,  pti, ptj), pose_i);
    dJJt_d_pose_j_num = numeric_jacobian_matrix_pose(@(pose_j) JJt(pose_i, pose_j, d,  pti, ptj), pose_j);
%     
   %  fprintf('dJJt_d_pose_i_num: \n');
   %  dJJt_d_pose_i_num
     fprintf('dJJt_d_pose_j_num: \n'); 
     dJJt_d_pose_j_num
    
    %% analyze jacobian 
    A = J*Jt; 
    K = size(J,2);
    y = reshape(A', [], 1); 
    T = size(y,1); 
    n = 6; % size(x,1); 
    dJ_dpi = zeros(T, n); 
    dJ_dpj = zeros(T, n); 
    for t=1:T %% row wise 
        i = int32(floor((t-1)/2)) + 1;
        j = t - (i-1)*2 ;
        for k = 1:K
            dJ_dpi(t,:) = dJ_dpi(t,:) + Jt(k,j)*dJ_dpose_i(K*(i-1)+k, :) + J(i,k)*dJt_dpose_i(j+(k-1)*2, :);  
            dJ_dpj(t,:) = dJ_dpj(t,:) + Jt(k,j)*dJ_dpose_j(K*(i-1)+k, :) + J(i,k)*dJt_dpose_j(j+(k-1)*2, :); 
        end
    end
    
    % fprintf('dJ_dpj: \n'); 
    % dJ_dpj
   
end

function A = JJt(pose_i, pose_j, d, pti, ptj)
    J = J_pose(pose_i, pose_j, d, pti, ptj);
    Jt = J'; 
    A = J*Jt; 
end
