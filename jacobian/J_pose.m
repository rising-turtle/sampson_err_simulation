
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