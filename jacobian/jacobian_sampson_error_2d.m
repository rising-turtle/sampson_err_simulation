%%
%   May 1, 2020, He Zhang, hzhang8@vcu.edu 
%   compute jacobian of the sampson error for 2d case 
%       se = ptj'*E*pti/sqrt(a^2 + b^2 + c^2 + d^2)
%       sse = J'*se, where J = [c, d, a, b]
%       a = (E*pti)_1, b = (E*pti)_2, c = (E'*ptj)_1, d = (E'*ptj)_2
%                           

function [dsse_dpi, dsse_dpj] = jacobian_sampson_error_2d(pose_i, pose_j, pti, ptj)

if nargin <= 0
    clc;
    pose_i = [1 2 3 sqrt(1-0.1^2-0.3^2-0.5^2) 0.1 -0.3 0.5]'; 
    pose_j = [0.1 0.2 0.3 sqrt(1-0.4^2-0.1^2-0.2^2) 0.4 0.1 -0.2]';
    pti = [0.1, 0.2]'; 
    ptj = [-0.3, 1.]';
    
end 

    %% numeric jacobian check 
    dse_dpose_i_num = numeric_jacobian_pose(@(pose_i) sampson_error_2d(pose_i, pose_j, pti, ptj), pose_i);   
    dse_dpose_i_num
    dse_dpose_j_num = numeric_jacobian_pose(@(pose_j) sampson_error_2d(pose_i, pose_j, pti, ptj), pose_j);   
    dse_dpose_j_num

    %% Jacobian of E w.r.t pose_i, pose_j 
    [dE_dpose_i, dE_dpose_j] = jacobian_essential_matrix(pose_i, pose_j, pti, ptj); 
    
    %% Essential matrix 
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
    
    %% compute 
    npti = [pti; 1]; 
    nptj = [ptj; 1]; 
    M = nptj*npti';
    lM = reshape(M', [], 1); 
    li = E*npti;
    lj = E'*nptj;
    a = li(1); b = li(2); c = lj(1); d = lj(2);
    
    %% se = Ue/Le 
    Ue = nptj'*E*npti;
    % Le = sqrt(a^2 + b^2 + c^2 + d^2); 
    Le = (a^2 + b^2 + c^2 + d^2); 
    
    %% dse_dUe
    dse_dpi = zeros(1,6); 
    dse_dpj = zeros(1,6); 
    
    for i=1:9
        dse_dpi = dse_dpi + lM(i)*dE_dpose_i(i, :);
        dse_dpj = dse_dpj + lM(i)*dE_dpose_j(i, :);  
    end
    dse_dpi = dse_dpi/Le; 
    dse_dpj = dse_dpj/Le;
    
    %% dse_dLe
    da_dpose_i = pti(1)*dE_dpose_i(1,:) + pti(2)*dE_dpose_i(2,:) + dE_dpose_i(3,:);
    db_dpose_i = pti(1)*dE_dpose_i(4,:) + pti(2)*dE_dpose_i(5,:) + dE_dpose_i(6,:);
    dc_dpose_i = ptj(1)*dE_dpose_i(1,:) + ptj(2)*dE_dpose_i(4,:) + dE_dpose_i(7,:);
    dd_dpose_i = ptj(1)*dE_dpose_i(2,:) + ptj(2)*dE_dpose_i(5,:) + dE_dpose_i(8,:);
    dle_dpose_i = -2*(a*da_dpose_i + b*db_dpose_i + c*dc_dpose_i + d*dd_dpose_i)/(Le^2);
    dse_dpi = dse_dpi + Ue*dle_dpose_i;
    % fprintf('dse_dpi: \n')
    % dse_dpi
    
    da_dpose_j = pti(1)*dE_dpose_j(1,:) + pti(2)*dE_dpose_j(2,:) + dE_dpose_j(3,:);
    db_dpose_j = pti(1)*dE_dpose_j(4,:) + pti(2)*dE_dpose_j(5,:) + dE_dpose_j(6,:);
    dc_dpose_j = ptj(1)*dE_dpose_j(1,:) + ptj(2)*dE_dpose_j(4,:) + dE_dpose_j(7,:);
    dd_dpose_j = ptj(1)*dE_dpose_j(2,:) + ptj(2)*dE_dpose_j(5,:) + dE_dpose_j(8,:);
    dle_dpose_j = -2*(a*da_dpose_j + b*db_dpose_j + c*dc_dpose_j + d*dd_dpose_j)/(Le^2);
    dse_dpj = dse_dpj + Ue*dle_dpose_j;
    % fprintf('dse_dpj: \n')
    % dse_dpj
    
    %% J 
    Jt = [c d a b]';
    se = Ue/Le; 
    dsse_dpi = zeros(4, 6);
    dsse_dpj = zeros(4, 6); 
    dsse_dpi = se*[dc_dpose_i; dd_dpose_i; da_dpose_i; db_dpose_i];
    dsse_dpj = se*[dc_dpose_j; dd_dpose_j; da_dpose_j; db_dpose_j];
    
    dsse_dpi = dsse_dpi + [c*dse_dpi; d*dse_dpi; a*dse_dpi; b*dse_dpi]; 
    dsse_dpj = dsse_dpj + [c*dse_dpj; d*dse_dpj; a*dse_dpj; b*dse_dpj]; 
    
    dsse_dpi
    dsse_dpj
    
end

function sse = sampson_error_2d(pose_i, pose_j, pti, ptj)
    
    %% Essential matrix 
    ti = pose_i(1:3); 
    qi = pose_i(4:7); % qw qx qy qz 
    tj = pose_j(1:3); 
    qj = pose_j(4:7);
    
    Ri = q2R(qi); 
    Rj = q2R(qj); 
    
    Rji = Rj'*Ri; 
    tij = Ri'*(tj-ti);
    Tij = skew(tij); 
    E = Rji*Tij;
    
    %% compute 
    npti = [pti; 1]; 
    nptj = [ptj; 1]; 

    li = E*npti;
    lj = E'*nptj;
    a = li(1); b = li(2); c = lj(1); d = lj(2);
    
    J = [c d a b]; 
    
    %% se = Ue/Le 
    Ue = nptj'*E*npti;
    % Le = sqrt(a^2 + b^2 + c^2 + d^2); 
    Le = a^2 + b^2 + c^2 + d^2;
    se = Ue/Le; 
    sse = J'*se;
end




