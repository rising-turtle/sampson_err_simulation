%%
%   April, 30, 2020, He Zhang, hzhang8@vcu.edu 
% compute jacobian of the matrix inv(JJt) 
%       where J =[de_dpti; de_dptj] 
%              de_dpti = [1 0 -ptj.y; -1 0 ptj.x] * Rj'* Ri * [d 0; 0 d; 0 0]
%              de_dptj = [0 -zj; zj 0]
%              

function [dJ_dpi, dJ_dpj] = jacobian_inv_JJt(pose_i, pose_j, d, pti, ptj)

if nargin <= 0
    pose_i = [1 2 3 sqrt(1-0.1^2-0.3^2-0.5^2) 0.1 -0.3 0.5]'; 
    pose_j = [0.1 0.2 0.3 sqrt(1-0.4^2-0.1^2-0.2^2) 0.4 0.1 -0.2]';
    d = 3;
    pti = [0.1, 0.2]'; 
    ptj = [-0.3, 1.]';
end
    J = J_pose(pose_i, pose_j, d, pti, ptj);
    Jt = J';    
    [dJJt_dpose_i, dJJt_dpose_j] = jacobian_JJt(pose_i, pose_j, d, pti, ptj);
    
    %% numeric result 
    diJJt_d_pose_i_num = numeric_jacobian_matrix_pose(@(pose_i) inv_JJt(pose_i, pose_j, d,  pti, ptj), pose_i);
    diJJt_d_pose_j_num = numeric_jacobian_matrix_pose(@(pose_j) inv_JJt(pose_i, pose_j, d,  pti, ptj), pose_j);
%     
     % fprintf('dJJt_d_pose_i_num: \n');
     % diJJt_d_pose_i_num
     % fprintf('diJJt_d_pose_j_num: \n'); 
     % diJJt_d_pose_j_num
    
    %% analyze jacobian 
    iJJt = inv_JJt(pose_i, pose_j, d, pti, ptj); 
    [m, n] = size(iJJt); 
    
    dJ_dpi = zeros(size(diJJt_d_pose_i_num));
    dJ_dpj = zeros(size(diJJt_d_pose_i_num));
    for c = 1:6
        dF_dxc = dJJt_dpose_i(:,c); 
        dF_dxc = reshape(dF_dxc', [m,n]); 
        dy_dx_m = -iJJt * dF_dxc * iJJt; 
        dJ_dpi(:,c) = reshape(dy_dx_m', [], 1); 
        dF_dxc = dJJt_dpose_j(:,c); 
        dF_dxc = reshape(dF_dxc', [m,n]); 
        dy_dx_m = -iJJt * dF_dxc * iJJt; 
        dJ_dpj(:,c) = reshape(dy_dx_m', [], 1); 
    end
    fprintf('dJ_dpj: \n');
    dJ_dpj
    
    %% test one matrix computation
    f11 = iJJt(1,1); f12 = iJJt(1,2); f21 = iJJt(2,1); f22 = iJJt(2,2); 
    F = -[f11^2,   f11*f21, f11*f12, f12*f21; 
         f11*f12, f11*f22, f12^2,   f12*f22; 
         f11*f21, f21^2,   f11*f22, f21*f22; 
         f12*f21, f21*f22, f12*f22, f22^2]; 
     dJ_dpj2 = F * dJJt_dpose_j; 
     fprintf('dJ_dpj2: \n'); 
     dJ_dpj2
    
end

function A = inv_JJt(pose_i, pose_j, d, pti, ptj)
    J = J_pose(pose_i, pose_j, d, pti, ptj);
    Jt = J'; 
    A = inv(J*Jt); 
end
