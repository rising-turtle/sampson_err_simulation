%%
%      April 29, 2020, He Zhang, hzhang8@vcu.edu
%       test rotation jacobian w.r.t quaternion  y = ARx, 
%           R is rotation matrix, q is R's quaternion 
%               y = Rp => dy_dtheta = -R*skew(p) 
%               y = R'p => dy_dtheta = skew(R'p)
%


function [dy_dx] = rotation_quaternion()

    x = [30, 0, 45]';  % euler engle 
    x = x*pi/180; 
    
    [R, dR] = e2R(x); 
    q = R2q(R); 
    p = [2, 0, 0]'; 
    % f = @(x) q2R(x)*p; 
    % f = @(x) q2R(x)'*p; 
    
    A = [0 1 -0.1; 
        -1 0 0.2];
    
    f = @(x) A*q2R(x)*p;
    
    fprintf('dy_dx_num: \n');
    [dy_dq_num] = numeric_jacobian_quaternion(f, q)
    
    %% analytical jacobian
    dy_dq = zeros(3, 3);
    
    % Sp = -skew_symmetric33(p);
    % dy_dq = R*Sp; 
    % pp = R' * p;
    % dy_dq = skew_symmetric33(pp);
    
    dy_dq = -A*R*skew_symmetric33(p);
    
    fprintf('dy_dx: \n'); 
    dy_dq
end
