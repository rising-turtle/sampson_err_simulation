%%
%      April 29, 2020, He Zhang, hzhang8@vcu.edu
%       test rotation jacobian w.r.t quaternion  y = ARx, 
%           R is rotation matrix, q is R's quaternion 
%               (dy_dq)_i = g'*dFi_dx + Fi.*dg_dx, 
%                   where (dy_dx)_i is the ith row of matrix dy_dx
%                   Fi is the ith row of F(x). 
%


function [dy_dx] = rotation_quaternion()

    x = [30, 0, 45]';  % euler engle 
    x = x*pi/180; 
    
    [R, dR] = e2R(x); 
    
    fprintf('dy_dx_num: \n');
    [dy_dx_num] = numeric_jacobian(@Ax, x)
    
    y = Ax(x); 
    
    m = size(y,1); 
    n = size(x,1); 
    
    dg_dx = 2*x; 
    g = gx(x);
    dy_dx = zeros(m, n); 
    for i=1:m
        dy_dx(i,:) = g'*dR(n*(i-1)+1:n*(i-1)+n, :) + R(i,:).*dg_dx';
    end
    
    fprintf('dy_dx: \n'); 
    dy_dx
end

function x = gx(x)
    x = x.^2; 
end

function y = Ax(x)

    A = e2R(x); 
    x = gx(x); 
    
    y = A*x; 

end