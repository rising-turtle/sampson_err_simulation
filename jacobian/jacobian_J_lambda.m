%%
%   March, 1, 2021, He Zhang, fuyinzh@gmail.com 
% compute jacobian of the matrix J w.r.t lambda 
%       where J =[de_dpti; de_dptj] 
%              de_dpti = [1 0 -ptj.y; -1 0 ptj.x] * Rj'* Ri * [d 0; 0 d; 0 0]
%              de_dptj = [0 -zj; zj 0]
%              d = 1/lambda 

function dJ_lambda = jacobian_J_lambda(pose_i, pose_j, d, pti, ptj)

if nargin <= 0
    
    pose_i = [1 2 3 sqrt(1-0.1^2-0.3^2-0.5^2) 0.1 -0.3 0.5]'; 
    pose_j = [0.1 0.2 0.3 sqrt(1-0.4^2-0.1^2-0.2^2) 0.4 0.1 -0.2]';
    d = 3;
    pti = [0.1, 0.2]'; 
    ptj = [-0.3, 1.]';
    
end
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
    e = [pj(2) - ptj(2)*pj(3); -pj(1) + ptj(1)*pj(3)]; 
    
    A = [0 1 -ptj(2); -1 0 ptj(1)]; 
    d1 = [d 0 0]'; 
    d2 = [0 d 0]'; 
    
    de_dpti = A*Rj'*Ri*[d1 d2]; 
    de_dptj = [0 -pj(3); pj(3) 0]; 
    
    %% jacobian for de_dpti and de_dptj     
    dd_lambda = -d^2; 
    lambda = 1/d; 
    r3 = Rji(3, :); 
    dzj_lambda = r3*[pti;1]*dd_lambda; % dpj(3)_lambda 
    dJ12 = A*Rj'*Ri*[1 0; 0 1; 0, 0]*dd_lambda; % first two cols 
    dJ34 = [0, -dzj_lambda; dzj_lambda, 0]; % second two cols 
    
    % check using numeric function 
    f_lambda = @(lambda) A*Rj'*Ri*[1, 0; 0, 1; 0, 0].*(1/lambda);
    [dJ_dlambda_12] = numeric_jacobian_matrix_scalar(f_lambda, lambda);

    f_lambda = @(lambda) [0 -1; 1 0].*([0, 0, 1]*Rj'*Ri*[pti; 1].*(1/lambda));
    [dJ_dlambda_34] = numeric_jacobian_matrix_scalar(f_lambda, lambda);
    
 if nargin <= 0   
    fprintf('dJ_dlambda_12_num \n');
    dJ_dlambda_12
    fprintf('dJ_dlambda_34_num \n');
    dJ_dlambda_34
    
    fprintf('dJ_dlambda_12 \n');
    dJ12  % = reshape(dJ12', [], 1) ;
    fprintf('dJ_dlambda_34 \n');
    dJ34 % = reshape(dJ34', [], 1) ;
 end
    dJ_lambda = [dJ12(1, :), dJ34(1,:), dJ12(2,:), dJ34(2,:)]';
    
    
    
end



