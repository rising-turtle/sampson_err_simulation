%%
%   May 1, 2020, He Zhang, hzhang8@vcu.edu 
%   compute jacobian of the essential matrix  
%       E = Rji*skew[tji]
%
%              

function [dE_dpose_i, dE_dpose_j] = jacobian_essential_matrix(pose_i, pose_j, pti, ptj)

if nargin <= 0
    clc;
    pose_i = [1 2 3 sqrt(1-0.1^2-0.3^2-0.5^2) 0.1 -0.3 0.5]'; 
    pose_j = [0.1 0.2 0.3 sqrt(1-0.4^2-0.1^2-0.2^2) 0.4 0.1 -0.2]';
    pti = [0.1, 0.2]'; 
    ptj = [-0.3, 1.]';
    
end 
    
    %% numeric jacobian 
    % dE_dpose_i_num = numeric_jacobian_matrix_pose(@(pose_i) essential_matrix(pose_i, pose_j), pose_i);   
    % dE_dpose_i_num
    % dE_dpose_j_num = numeric_jacobian_matrix_pose(@(pose_j) essential_matrix(pose_i, pose_j), pose_j);   
    % dE_dpose_j_num
    
    %% analytical jacobian 
    
    ti = pose_i(1:3); 
    qi = pose_i(4:7); % qw qx qy qz 
    tj = pose_j(1:3); 
    qj = pose_j(4:7);
    
    Ri = q2R(qi); 
    Rj = q2R(qj); 
    
    Rji = Rj'*Ri; 
    tij = Ri'*(tj-ti);
    % tji = Rj'*(ti-tj);
    Tij = skew(tij); 
    E = Rji*Tij; 
    dEc1_dtheta_i = -Rj'*Ri*skew_symmetric33(Tij(:,1)); 
    dEc2_dtheta_i = -Rj'*Ri*skew_symmetric33(Tij(:,2));
    dEc3_dtheta_i = -Rj'*Ri*skew_symmetric33(Tij(:,3));
    dEc1_dtheta_j = skew_symmetric33(E(:,1));
    dEc2_dtheta_j = skew_symmetric33(E(:,2));
    dEc3_dtheta_j = skew_symmetric33(E(:,3));
    
    %% dE_dtij, tij = Ri'*(tj-ti) 
    dE_dtij = jacobian_skew_t(Rji, tij); 
    
    %% dtij_dpose_i, dtij_dpose_j 
    dtij_dtheta_i = skew_symmetric33(tij); 
    dtij_dtheta_j = zeros(3,3); 
    dtij_dti = -Ri'; 
    dtij_dtj = Ri'; 
    dtij_dpose_i = [dtij_dti dtij_dtheta_i]; 
    dtij_dpose_j = [dtij_dtj dtij_dtheta_j];
    dE_dpose_i = dE_dtij*dtij_dpose_i; 
    dE_dpose_j = dE_dtij*dtij_dpose_j;    
    
    %% add derivative with regard to Rji part 
    dE_dpose_i(1,4:6) = dE_dpose_i(1,4:6) + dEc1_dtheta_i(1,:);
    dE_dpose_i(2,4:6) = dE_dpose_i(2,4:6) + dEc2_dtheta_i(1,:);
    dE_dpose_i(3,4:6) = dE_dpose_i(3,4:6) + dEc3_dtheta_i(1,:);
    dE_dpose_i(4,4:6) = dE_dpose_i(4,4:6) + dEc1_dtheta_i(2,:);
    dE_dpose_i(5,4:6) = dE_dpose_i(5,4:6) + dEc2_dtheta_i(2,:);
    dE_dpose_i(6,4:6) = dE_dpose_i(6,4:6) + dEc3_dtheta_i(2,:);
    dE_dpose_i(7,4:6) = dE_dpose_i(7,4:6) + dEc1_dtheta_i(3,:);
    dE_dpose_i(8,4:6) = dE_dpose_i(8,4:6) + dEc2_dtheta_i(3,:);
    dE_dpose_i(9,4:6) = dE_dpose_i(9,4:6) + dEc3_dtheta_i(3,:);
    % dE_dpose_i
    
    dE_dpose_j(1,4:6) = dE_dpose_j(1,4:6) + dEc1_dtheta_j(1,:);
    dE_dpose_j(2,4:6) = dE_dpose_j(2,4:6) + dEc2_dtheta_j(1,:);
    dE_dpose_j(3,4:6) = dE_dpose_j(3,4:6) + dEc3_dtheta_j(1,:);
    dE_dpose_j(4,4:6) = dE_dpose_j(4,4:6) + dEc1_dtheta_j(2,:);
    dE_dpose_j(5,4:6) = dE_dpose_j(5,4:6) + dEc2_dtheta_j(2,:);
    dE_dpose_j(6,4:6) = dE_dpose_j(6,4:6) + dEc3_dtheta_j(2,:);
    dE_dpose_j(7,4:6) = dE_dpose_j(7,4:6) + dEc1_dtheta_j(3,:);
    dE_dpose_j(8,4:6) = dE_dpose_j(8,4:6) + dEc2_dtheta_j(3,:);
    dE_dpose_j(9,4:6) = dE_dpose_j(9,4:6) + dEc3_dtheta_j(3,:);
    % dE_dpose_j
end

function dE_dt = jacobian_skew_t(R, t)

    dE_dt = zeros(9, 3); 
    
    dE_dtx = R*[0 0 0; 0 0 -1; 0 1 0]; 
    dE_dty = R*[0 0 1; 0 0 0; -1 0 0];
    dE_dtz = R*[0 -1 0; 1 0 0; 0 0 0]; 
    dE_dt(:,1) = reshape(dE_dtx', [], 1); 
    dE_dt(:,2) = reshape(dE_dty', [], 1); 
    dE_dt(:,3) = reshape(dE_dtz', [], 1); 

end

function E = essential_matrix(pose_i, pose_j)
    ti = pose_i(1:3); 
    qi = pose_i(4:7); % qw qx qy qz 
    tj = pose_j(1:3); 
    qj = pose_j(4:7);
    
    Ri = q2R(qi); 
    Rj = q2R(qj); 
    
    Rji = Rj'*Ri; 
    % tji = Rj'*(ti-tj);
    tij = Ri'*(tj-ti);
    
    E = Rji*skew(tij); 
end


function E = essential_matrix_tmp(pose_i, pose_j, Rji)
    ti = pose_i(1:3); 
    qi = pose_i(4:7); % qw qx qy qz 
    tj = pose_j(1:3); 
    qj = pose_j(4:7);
    
    Ri = q2R(qi); 
    % Rj = q2R(qj); 
    
    % Rji = Rj'*Ri; 
    % tji = Rj'*(ti-tj);
    tij = Ri'*(tj-ti);
    
    E = Rji*skew(tij); 
end