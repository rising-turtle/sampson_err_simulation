%%
%    April 29, 2020, He Zhang, hzhang8@vcu.edu
%       numeric jacobians where y = Rx is a matrix
%           output dy_dq
%

function [dy_dtheta] = numeric_jacobian_quaternion(f,q)

a = [1, 2, 3]'; 
if nargin <= 0   
    f = @(x)q2R(x)*a; 
    x = [1., 2.,0.5]';
    R = e2R(x); 
    q = R2q(R); 
end

R = q2R(q); 
de = 1e-6; 
y = f(q);
n = 3;
dy_dtheta = zeros(size(), 3); 

for i = 1:3
    z = zeros(n,1); 
    z(i) = de; 
    dq = delta_q(z);
    x2 = qProd(q, dq);
    y2 = f(x2);
    
    dq = delta_q(-z); 
    x1 = qProd(q, dq); 
    y1 = f(x1); 
    dy_dtheta(:,i) = (y2-y1)/(2*de); 
end

%% test the result 
% dtheta = 1e-3 * ones(3,1); 
% new_q = qProd(q, delta_q(dtheta));
% new_y = f(new_q); 
% 
% % estimated 
% new_y_est = f(q) + dy_dtheta * dtheta; 
% fprintf('new_y: %f new_y_est: %f \r\n', new_y, new_y_est); 

end





