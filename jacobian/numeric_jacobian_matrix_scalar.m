%%
%    March 1, 2021, He Zhang, fuyinzh@vcu.edu
%       numeric jacobians where A = f(x) is a matrix, 
%           x is a scalar vector
%           y = reshape(A) row wise
%           output dy_dx
%

function [dy_dx] = numeric_jacobian_matrix_scalar(f,x)

if nargin <= 0
    % f = @(x)e2R(x); 
    % x = [1., 2.,0.5]';
    % [R, dR] = e2R(x);
end

de = 1e-6; 

n = size(x,1);
A = f(x); 
y = reshape(A', [], 1); 
m = size(y,1); 
dy_dx = zeros(m, n); 

for i = 1:n
    z = zeros(n,1); 
    z(i) = de;
    x1 = x - z; 
    x2 = x + z; 
    A2 = f(x2);
    y2 = reshape(A2', [], 1);
    A1 = f(x1); 
    y1 = reshape(A1', [], 1);  
    dy_dx(:,i) = (y2-y1)/(2*de); 
end



end

