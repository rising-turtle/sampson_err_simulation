%%
%   March, 1, 2021, He Zhang, fuyinzh@gmail.com  
% compute jacobian of the matrix JJt 
%       where J =[de_dpti; de_dptj] 
%              de_dpti = [1 0 -ptj.y; -1 0 ptj.x] * Rj'* Ri * [d 0; 0 d; 0 0]
%              de_dptj = [0 -zj; zj 0]
%              

function dJJt_lambda = jacobian_JJt(pose_i, pose_j, d, pti, ptj)

if nargin <= 0
    pose_i = [1 2 3 sqrt(1-0.1^2-0.3^2-0.5^2) 0.1 -0.3 0.5]'; 
    pose_j = [0.1 0.2 0.3 sqrt(1-0.4^2-0.1^2-0.2^2) 0.4 0.1 -0.2]';
    d = 3;
    pti = [0.1, 0.2]'; 
    ptj = [-0.3, 1.]';
end
    J = J_pose(pose_i, pose_j, d, pti, ptj);
    Jt = J'; 
    dJ_lambda = jacobian_J_lambda(pose_i, pose_j, d, pti, ptj);
    dJt_lambda = jacobian_Jt_lambda(pose_i, pose_j, d, pti, ptj);
    
    lambda = 1/d;
    %% numeric result 
    dJJt_d_lambda_num = numeric_jacobian_matrix_scalar(@(lambda) JJt(pose_i, pose_j, lambda,  pti, ptj), lambda);
    
   %  fprintf('dJJt_d_pose_i_num: \n');
   %  dJJt_d_pose_i_num
    
    %% analyze jacobian 
    A = J*Jt; 
    K = size(J,2);
    y = reshape(A', [], 1); 
    T = size(y,1); 
    n = 1; % size(x,1); 
    dJJt_lambda = zeros(T, 1); 
    for t=1:T %% row wise 
        i = int32(floor((t-1)/2)) + 1;
        j = t - (i-1)*2 ;
        for k = 1:K
            dJJt_lambda(t,:) = dJJt_lambda(t,:) + Jt(k,j)*dJ_lambda(K*(i-1)+k, :) + J(i,k)*dJt_lambda(j+(k-1)*2, :);  
        end
    end
if nargin <= 0
    fprintf('dJJt_d_pose_j_num: \n'); 
    dJJt_d_lambda_num
    fprintf('dJJt_lambda: \n'); 
    dJJt_lambda
end

end

function A = JJt(pose_i, pose_j, lambda, pti, ptj)
    J = J_pose(pose_i, pose_j, 1/lambda, pti, ptj);
    Jt = J'; 
    A = J*Jt; 
end
