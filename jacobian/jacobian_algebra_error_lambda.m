%%
%   March 1, 2021, He Zhang, fuyinzh@gmail.com
%       jacobian of algebra error 
%           de_d_lambda
%

function de_d_lambda= jacobian_algebra_error_lambda(pose_i, pose_j, d, pti, ptj)

if nargin <= 0
    pose_i = [1 2 3 sqrt(1-0.1^2-0.3^2-0.5^2) 0.1 -0.3 0.5]'; 
    pose_j = [0.1 0.2 0.3 sqrt(1-0.4^2-0.1^2-0.2^2) 0.4 0.1 -0.2]';
    d = 3;
    pti = [0.1, 0.2]'; 
    ptj = [-0.3, 1.]';
end

%% numeric jacobian of algebra error 
lambda = 1/d; 
de_d_lambda_num = numeric_jacobian_matrix_scalar(@(lambda) algebra_error(pose_i, pose_j, 1/lambda, pti, ptj), lambda); 

if nargin <= 0
    fprintf('de_d_lambda_num: \n');
    de_d_lambda_num
end

%% analytical
ti = pose_i(1:3);
qi = pose_i(4:7); % qw qx qy qz
tj = pose_j(1:3);
qj = pose_j(4:7);

Ri = q2R(qi);
Rj = q2R(qj);
Rji = Rj'*Ri; 
r1 = Rji(1, :); 
r2 = Rji(2, :);
r3 = Rji(3, :); 
pj = Rj'*Ri*[pti; 1]*d + Rj'*(ti-tj);
e = algebra_error(pose_i, pose_j, d, pti, ptj); % e = [ y - v*z; -x + u * z]

dd_d_lambda = -d^2; 
dx_d_lambda = r1 * [pti; 1] * dd_d_lambda;  
dy_d_lambda = r2 * [pti; 1] * dd_d_lambda;   
dz_d_lambda = r3 * [pti; 1] * dd_d_lambda;   

de1_d_lambda = dy_d_lambda - ptj(2) * dz_d_lambda; 
de2_d_lambda = - dx_d_lambda + ptj(1) * dz_d_lambda; 
de_d_lambda = [de1_d_lambda; de2_d_lambda]; 
if nargin <= 0
    fprintf('de_d_lambda: \n');
    de_d_lambda
end


end

