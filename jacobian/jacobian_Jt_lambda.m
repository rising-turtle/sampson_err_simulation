%%
%   March, 1, 2021, He Zhang, fuyinzh@gmail.com 
% compute jacobian of the matrix Jt w.r.t lambda 
%       where J =[de_dpti; de_dptj] 
%              de_dpti = [1 0 -ptj.y; -1 0 ptj.x] * Rj'* Ri * [d 0; 0 d; 0 0]
%              de_dptj = [0 -zj; zj 0]
%              d = 1/lambda 

function dJt_lambda = jacobian_Jt_lambda(pose_i, pose_j, d, pti, ptj)

if nargin <= 0
    
    pose_i = [1 2 3 sqrt(1-0.1^2-0.3^2-0.5^2) 0.1 -0.3 0.5]'; 
    pose_j = [0.1 0.2 0.3 sqrt(1-0.4^2-0.1^2-0.2^2) 0.4 0.1 -0.2]';
    d = 3;
    pti = [0.1, 0.2]'; 
    ptj = [-0.3, 1.]';
    
end
    J = jacobian_J_lambda(pose_i, pose_j, d, pti, ptj); 
    dJt_lambda = [J(1), J(5), J(2), J(6), J(3), J(7), J(4), J(8)]';
    
end



