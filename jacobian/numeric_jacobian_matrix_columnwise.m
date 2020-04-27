%%
%    April 26, 2020, He Zhang, hzhang8@vcu.edu
%       numeric jacobians where A = f(x) is a matrix
%           y = reshape(A) column wise
%           output dy_dx
%

function [dy_dx] = numeric_jacobian_matrix_columnwise(f,x)

if nargin <= 0
    f = @(x)e2R(x); 
    x = [1., 2.,0.5]';
    [R, dR] = e2R(x);
end

de = 1e-6; 

n = size(x,1);
A = f(x); 
y = reshape(A, [], 1); 
m = size(y,1); 
dy_dx = zeros(m, n); 

for i = 1:size(x,1)
    z = zeros(n,1); 
    z(i) = de; 
    x2 = x + z;
    A2 = f(x2);
    y2 = reshape(A2, [], 1); 
    x1 = x - z; 
    A1 = f(x1); 
    y1 = reshape(A1, [], 1);  
    dy_dx(:,i) = (y2-y1)/(2*de); 
    
end



end

