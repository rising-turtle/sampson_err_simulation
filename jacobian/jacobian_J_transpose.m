%%
%   April, 30, 2020, He Zhang, hzhang8@vcu.edu 
% compute jacobian of the matrix J' 
%       where J =[de_dpti; de_dptj] 
%              de_dpti = [1 0 -ptj.y; -1 0 ptj.x] * Rj'* Ri * [d 0; 0 d; 0 0]
%              de_dptj = [0 -zj; zj 0]
%              

function [dJ_dpi, dJ_dpj] = jacobian_J_transpose(pose_i, pose_j, d, pti, ptj)

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
    e = [pj(2) - ptj(2)*pj*(3); -pj(1) + ptj(1)*pj(3)]; 
    
    A = [0 1 -ptj(2); -1 0 ptj(1)]; 
    d1 = [d 0 0]'; 
    d2 = [0 d 0]'; 
    
    de_dpti = A*Rj'*Ri*[d1 d2]; 
    de_dptj = [0 -pj(3); pj(3) 0]; 
    
    %% jacobian for de_dpti
    J1 = de_dpti(:,1); 
    J2 = de_dpti(:,2); 
    dJ1_dtheta_i = -A*Rj'*Ri*skew_symmetric33(d1);
    dJ1_dtheta_j = A*skew_symmetric33(Rj'*Ri*d1); 
    dJ1_dti = zeros(2,3); 
    dJ1_dtj = zeros(2,3);
    
    dJ2_dtheta_i = -A*Rj'*Ri*skew_symmetric33(d2);
    dJ2_dtheta_j = A*skew_symmetric33(Rj'*Ri*d2); 
    dJ2_dti = zeros(2,3); 
    dJ2_dtj = zeros(2,3);
    
    %% jacobian for de_dptj 
    % Rji = Rj'*Ri; 
    % tji = Rj'*(ti-tj);
    % pj = Rj'*Ri*[pti; 1]*d + tji; 
    dpj_dtheta_i = -Rj'*Ri*skew_symmetric33([pti; 1]*d);
    dpj_dti = Rj';
    dpj_dtheta_j = skew_symmetric33(Rj'*Ri*[pti; 1]*d) + skew_symmetric33(tji);
    dpj_dtj = -Rj';
 
    r3 = [0 0 1]; 
    dzj_dtheta_i = r3*dpj_dtheta_i; 
    dzj_dti = r3*dpj_dti; 
    
    dJ34_dpose_i = [[0 0 0 0 0 0]; 
                    [-dzj_dti -dzj_dtheta_i];
                    [dzj_dti  dzj_dtheta_i];
                    [0 0 0 0 0 0]];
    
    dzj_dtheta_j = r3*dpj_dtheta_j; 
    dzj_dtj = r3*dpj_dtj;
    dJ34_dpose_j = [[0 0 0 0 0 0]; 
                    [-dzj_dtj -dzj_dtheta_j];
                    [dzj_dtj  dzj_dtheta_j];
                    [0 0 0 0 0 0]];
    
    % check using numeric function
    
    % finally, combine the result 
    dJ_dpi = [[dJ1_dti(1,:) dJ1_dtheta_i(1,:)]; 
              [dJ1_dti(2,:) dJ1_dtheta_i(2,:)];
              [dJ2_dti(1,:) dJ2_dtheta_i(1,:)];
              [dJ2_dti(2,:) dJ2_dtheta_i(2,:)];
              dJ34_dpose_i(1,:);
              dJ34_dpose_i(3,:);
              dJ34_dpose_i(2,:);
              dJ34_dpose_i(4,:)
              ];
    dJ_dpj = [[dJ1_dtj(1,:) dJ1_dtheta_j(1,:)]; 
              [dJ1_dtj(2,:) dJ1_dtheta_j(2,:)];
              [dJ2_dtj(1,:) dJ2_dtheta_j(1,:)];
              [dJ2_dtj(2,:) dJ2_dtheta_j(2,:)];
              dJ34_dpose_j(1,:);
              dJ34_dpose_j(3,:);
              dJ34_dpose_j(2,:);
              dJ34_dpose_j(4,:)
              ];   
  % check using numeric function 
%     dJ_d_pose_i_num = numeric_jacobian_matrix_pose(@(pose_i) J_pose(pose_i, pose_j, d, pti, ptj), pose_i);
%     dJ_d_pose_j_num = numeric_jacobian_matrix_pose(@(pose_j) J_pose(pose_i, pose_j, d, pti, ptj), pose_j);
%     
%     fprintf('dJ_d_pose_i_num: \n');
%     dJ_d_pose_i_num
%     fprintf('dJ_d_pose_j_num: \n'); 
%     dJ_d_pose_j_num

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
    J = [de_dpti de_dptj]';

end
