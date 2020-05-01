%%
%   April 30, 2020, He Zhang, hzhang8@vcu.edu
%       jacobian of algebra error 
%           [de_dpose_i, de_dpose_j]
%

function [de_dpose_i, de_dpose_j] = jacobian_algebra_error(pose_i, pose_j, d, pti, ptj)

if nargin <= 0
    pose_i = [1 2 3 sqrt(1-0.1^2-0.3^2-0.5^2) 0.1 -0.3 0.5]'; 
    pose_j = [0.1 0.2 0.3 sqrt(1-0.4^2-0.1^2-0.2^2) 0.4 0.1 -0.2]';
    d = 3;
    pti = [0.1, 0.2]'; 
    ptj = [-0.3, 1.]';
end

%% numeric jacobian of algebra error 
de_dpose_i_num = numeric_jacobian_pose(@(pose_i) algebra_error(pose_i, pose_j, d, pti, ptj), pose_i); 
de_dpose_j_num = numeric_jacobian_pose(@(pose_j) algebra_error(pose_i, pose_j, d, pti, ptj), pose_j); 

% fprintf('de_dpose_i_num: \n');
% de_dpose_i_num
% fprintf('de_dpose_j_num: \n'); 
% de_dpose_j_num

%% analytical
ti = pose_i(1:3);
qi = pose_i(4:7); % qw qx qy qz
tj = pose_j(1:3);
qj = pose_j(4:7);

Ri = q2R(qi);
Rj = q2R(qj);
pj = Rj'*Ri*[pti; 1]*d + Rj'*(ti-tj);

e = algebra_error(pose_i, pose_j, d, pti, ptj);

de_dpj = -skew([ptj; 1]); 
dpj_dti = Rj';
dpj_dtheta_i = -Rj'*Ri*skew_symmetric33([pti; 1]*d);
dpj_dtj = -Rj'; 
dpj_dtheta_j = skew_symmetric33(Rj'*Ri*[pti; 1]*d) + skew_symmetric33(Rj'*(ti-tj));

ext_m = [1 0 0; 0 1 0];
de_dpose_i = ext_m * de_dpj * [dpj_dti dpj_dtheta_i];
de_dpose_j = ext_m * de_dpj * [dpj_dtj dpj_dtheta_j]; 
% fprintf('de_dpose_i: \n');
% de_dpose_i
% fprintf('de_dpose_j: \n'); 
% de_dpose_j


end

