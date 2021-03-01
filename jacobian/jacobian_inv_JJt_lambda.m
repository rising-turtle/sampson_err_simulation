%%
%   March, 1, 2021, He Zhang, fuyinzh@gmail.com
% compute jacobian of the matrix inv(JJt) 
%       where J =[de_dpti; de_dptj] 
%              de_dpti = [1 0 -ptj.y; -1 0 ptj.x] * Rj'* Ri * [d 0; 0 d; 0 0]
%              de_dptj = [0 -zj; zj 0]
%              

function diJJt_d_lambda_2 = jacobian_inv_JJt_lambda(pose_i, pose_j, d, pti, ptj)

if nargin <= 0
    pose_i = [1 2 3 sqrt(1-0.1^2-0.3^2-0.5^2) 0.1 -0.3 0.5]'; 
    pose_j = [0.1 0.2 0.3 sqrt(1-0.4^2-0.1^2-0.2^2) 0.4 0.1 -0.2]';
    d = 3;
    pti = [0.1, 0.2]'; 
    ptj = [-0.3, 1.]';
end
    J = J_pose(pose_i, pose_j, d, pti, ptj);
    Jt = J';    
    dJJt_d_lambda = jacobian_JJt_lambda(pose_i, pose_j, d, pti, ptj);
    
    lambda = 1/d;
    %% numeric result 
    diJJt_d_lambda_num = numeric_jacobian_matrix_scalar(@(lambda) inv_JJt(pose_i, pose_j, lambda,  pti, ptj), lambda);
    % fprintf('diJJt_d_lambda_num: \n');
    % diJJt_d_lambda_num
    
     % fprintf('dJJt_d_pose_i_num: \n');
     % diJJt_d_pose_i_num
     % fprintf('diJJt_d_pose_j_num: \n'); 
     % diJJt_d_pose_j_num
    
    %% analyze jacobian 
    iJJt = inv_JJt(pose_i, pose_j, lambda, pti, ptj); 
    [m, n] = size(iJJt); 
    
    diJJt_d_lambda = zeros(size(diJJt_d_lambda_num));
   
    dF_dxc = dJJt_d_lambda(:,1); 
    dF_dxc = reshape(dF_dxc', [m,n]); 
    dy_dx_m = -iJJt * dF_dxc * iJJt; 
    diJJt_d_lambda = reshape(dy_dx_m', [], 1); 
    % fprintf('diJJt_d_lambda: \n'); 
    % diJJt_d_lambda
    
    %% test one matrix computation
    f11 = iJJt(1,1); f12 = iJJt(1,2); f21 = iJJt(2,1); f22 = iJJt(2,2); 
    F = -[f11^2,   f11*f21, f11*f12, f12*f21; 
         f11*f12, f11*f22, f12^2,   f12*f22; 
         f11*f21, f21^2,   f11*f22, f21*f22; 
         f12*f21, f21*f22, f12*f22, f22^2]; 
     diJJt_d_lambda_2 = F * dJJt_d_lambda; 
     % fprintf('diJJt_d_lambda_2: \n'); 
     % diJJt_d_lambda_2
    
end

function A = inv_JJt(pose_i, pose_j, lambda, pti, ptj)
    J = J_pose(pose_i, pose_j, 1/lambda, pti, ptj);
    Jt = J'; 
    A = inv(J*Jt); 
end
