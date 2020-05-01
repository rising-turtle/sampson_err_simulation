%%
%   April, 30, 2020, He Zhang, hzhang8@vcu.edu 
%   compute jacobian of the sampson error  
%       se = -J'inv(JJ')e
%       where J =[de_dpti; de_dptj] 
%              de_dpti = [1 0 -ptj.y; -1 0 ptj.x] * Rj'* Ri * [d 0; 0 d; 0 0]
%              de_dptj = [0 -zj; zj 0]
%              e is algebra error, cross_product(ptj, pj) 
%                   where pj = Rji*[pti;1]*d + tji
%              

function [dse_dpi, dse_dpj] = jacobian_sampson_error(pose_i, pose_j, d, pti, ptj)

if nargin <= 0
    clc;
    pose_i = [1 2 3 sqrt(1-0.1^2-0.3^2-0.5^2) 0.1 -0.3 0.5]'; 
    pose_j = [0.1 0.2 0.3 sqrt(1-0.4^2-0.1^2-0.2^2) 0.4 0.1 -0.2]';
    d = 3;
    pti = [0.1, 0.2]'; 
    ptj = [-0.3, 1.]';
    
end 

    %% dz_dpi, z = inv(JJ')*e
    dz_dpi = zeros(2,6);
    dz_dpj = zeros(2,6); 
    [dA_dpose_i, dA_dpose_j] = jacobian_inv_JJt(pose_i, pose_j, d, pti, ptj);
    % de_dpose_i = numeric_jacobian_pose(@(pose_i) algebra_error(pose_i, pose_j, d, pti, ptj), pose_i); 
    % de_dpose_j = numeric_jacobian_pose(@(pose_j) algebra_error(pose_i, pose_j, d, pti, ptj), pose_j); 
    [de_dpose_i, de_dpose_j] = jacobian_algebra_error(pose_i, pose_j, d, pti, ptj);
    
    e = algebra_error(pose_i, pose_j, d, pti, ptj); 
    J = J_pose(pose_i, pose_j, d, pti, ptj);
    A = inv(J*J');
    k = size(A, 2);
    for i=1:2
        for j=1:size(e,1)
            dz_dpi(i,:) = dz_dpi(i,:) + e(j)*dA_dpose_i(k*(i-1)+j, :); 
            dz_dpj(i,:) = dz_dpj(i,:) + e(j)*dA_dpose_j(k*(i-1)+j, :); 
        end
        for j=1:size(A,2)
            dz_dpi(i,:) = dz_dpi(i,:) + A(i,j)*de_dpose_i(j, :);
            dz_dpj(i,:) = dz_dpj(i,:) + A(i,j)*de_dpose_j(j, :);
        end
    end
    % fprintf('dz_dpi: \n')
    % dz_dpi
    
    % check using numeric function 
    % dz_d_pose_i_num = numeric_jacobian_matrix_pose(@(pose_i) invJJt_e(pose_i, pose_j, d, pti, ptj), pose_i);
    % dJ_d_pose_j_num = numeric_jacobian_matrix_pose(@(pose_j) sampson_error(pose_i, pose_j, d, pti, ptj), pose_j);
    
    % fprintf('dz_d_pose_i_num: \n');
    % dz_d_pose_i_num
    % fprintf('dJ_d_pose_j_num: \n'); 
    
    % check using numeric function 
    dS_d_pose_i_num = numeric_jacobian_matrix_pose(@(pose_i) sampson_error(pose_i, pose_j, d, pti, ptj), pose_i);
    dS_d_pose_j_num = numeric_jacobian_matrix_pose(@(pose_j) sampson_error(pose_i, pose_j, d, pti, ptj), pose_j);
    
    fprintf('dS_d_pose_i_num: \n');
    dS_d_pose_i_num
    fprintf('dS_d_pose_j_num: \n'); 
    dS_d_pose_j_num
    %% initialize
    m = 4; 
    n = 6; 
    dse_dpi = zeros(4, 6); 
    dse_dpj = zeros(4, 6); 
    
    %% dJ'_dpi
    [dJt_dpose_i, dJt_dpose_j] = jacobian_J_transpose(pose_i, pose_j, d, pti, ptj);
    ze = A * e; 
    Jt = J_pose(pose_i, pose_j, d, pti, ptj)';
    
    k = size(Jt, 2);
    for i=1:m
        for j=1:size(ze,1)
            dse_dpi(i,:) = dse_dpi(i,:) + ze(j)*dJt_dpose_i(k*(i-1)+j, :);
            dse_dpj(i,:) = dse_dpj(i,:) + ze(j)*dJt_dpose_j(k*(i-1)+j, :);
        end
        for j=1:size(Jt,2)
            dse_dpi(i,:) = dse_dpi(i,:) + Jt(i,j)*dz_dpi(j, :);
            dse_dpj(i,:) = dse_dpj(i,:) + Jt(i,j)*dz_dpj(j, :);
        end
    end
    fprintf('dse_dpi: \n')
    dse_dpi
    fprintf('dse_dpj: \n')
    dse_dpj
end

function z = invJJt_e(pose_i, pose_j, d, pti, ptj)
     J = J_pose(pose_i, pose_j, d, pti, ptj); 
     A = inv(J*J'); 
     e = algebra_error(pose_i, pose_j, d, pti, ptj); 
     z = A*e;
end

function se = sampson_error(pose_i, pose_j, d, pti, ptj)
    
    J = J_pose(pose_i, pose_j, d, pti, ptj); 
    A = inv(J*J');
    e = algebra_error(pose_i, pose_j, d, pti, ptj); 
    ze = A*e; 
    % step by step 
    se = J'*ze; 

end

function J = J_pose(pose_i, pose_j, d, pti, ptj)
    ti = pose_i(1:3); 
    qi = pose_i(4:7); % qw qx qy qz 
    tj = pose_j(1:3); 
    qj = pose_j(4:7);
    
    Ri = q2R(qi); 
    Rj = q2R(qj); 
    
    Rji = Rj'*Ri; 
    tji = Rj'*(ti-tj);
    pj = Rji*[pti; 1]*d + tji; 
    
    %% define sampson error related terms     
    A = [0 1 -ptj(2); -1 0 ptj(1)]; 
    d1 = [d 0 0]'; 
    d2 = [0 d 0]'; 
    
    de_dpti = A*Rj'*Ri*[d1 d2]; 
    de_dptj = [0 -pj(3); pj(3) 0]; 
    J = [de_dpti de_dptj];

end


