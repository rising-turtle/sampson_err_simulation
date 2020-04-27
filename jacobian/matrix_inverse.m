%%
%      April 26, 2020, He Zhang, hzhang8@vcu.edu
%       test jacobians A = (F(x))-1, y = reshape(A', [], 1)
%           where F is a (mxm) semi positive definite matrix 
%               (dA_dx)ij = Fi.*dGj_dx + Gj.*dFi_dx, 
%                   where Fi is the ith row of matrix F
%                   Gj is the ith column of matrix G. 
%

function [dy_dx] = matrix_inverse()

    x = [30, 20, 45]';  % euler engle 
    x = x*pi/180; 

    [dy_dx_num] = numeric_jacobian_matrix(@Fv, x);
    F = Fx(x); 
    iF = inv(F);
    [m,n] = size(F); 
    fprintf('dy_dx_num: \n');
    % dy_dx_num = reshape(dy_dx_num', [m, n]);
    [dF_dx] = numeric_jacobian_matrix(@Fx, x); 
    
    dy_dx = zeros(size(dy_dx_num));
    
    for c = 1:n
        dF_dxc = dF_dx(:,c); 
        dF_dxc = reshape(dF_dxc', [m,n]); 
        dy_dx_m = -iF * dF_dxc * iF; 
        dy_dx(:,c) = reshape(dy_dx_m', [], 1); 
    end
    fprintf('dy_dx: \n');
    dy_dx
end

function G = Gx(x)
    G = x*x'; 
end

function F = Fx(x)
    n = size(x,1); 
    X = Gx(x) + eye(n);
    F = X*X'; 
end

function Fv = Fv(x)
    F = Fx(x); 
    Fv = inv(F);
end
