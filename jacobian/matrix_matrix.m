%%
%      April 26, 2020, He Zhang, hzhang8@vcu.edu
%       test jacobians A = F(x)*G(x), y = reshape(A', [], 1)
%           where F is a (mxk) matrix whereas G is a (kxn) vector
%                 A is a (mxn) matrix, y = reshape(A', [], 1)
%               (dA_dx)ij = Fi.*dGj_dx + Gj.*dFi_dx, 
%                   where Fi is the ith row of matrix F
%                   Gj is the ith column of matrix G. 
%


function [dy_dx] = matrix_matrix()

    x = [30, 0, 45]';  % euler engle 
    x = x*pi/180; 
    
    fprintf('dy_dx_num: \n');
    [dy_dx_num] = numeric_jacobian_matrix(@Ax, x);
    
    [dG_dx] = numeric_jacobian_matrix_columnwise(@Gx, x); % column wise
    [dF_dx] = numeric_jacobian_matrix(@e2R, x); % row wise
    
    F = e2R(x);
    G = Gx(x);
    m = size(F,1);
    k = size(F,2);
    A = Ax(x); 
    
    y = reshape(A', [], 1); 
    T = size(y,1); 
    n = size(x,1); 
    
    dy_dx = zeros(T, n); 
    for t=1:T %% row wise 
        i = int32(floor((t-1)/n)) + 1;
        j = t - (i-1)*n ;
        dy_dx(t,:) = G(j,:)*dF_dx(k*(i-1)+1:k*i, :) + F(i,:)*dG_dx(k*(j-1)+1:k*j, :);  
    end
    
    fprintf('dy_dx: \n'); 
    dy_dx
end

function G = Gx(x)
    G = x*x'; 
end

function A = Ax(x)

    F = e2R(x); 
    G = Gx(x); 
    
    A = F*G; 

end