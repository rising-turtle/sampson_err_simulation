%%
%    April 26, 2020, He Zhang, hzhang8@vcu.edu
%       numeric jacobians 
%           y = f(x)
%           output dy_dx
%

function [dy_dx] = numeric_jacobian(f,x)

if nargin <= 0
    f = @(x) 2*x'*x+3*norm(x); 
    x = [1., 2.,0.5]';
end

de = 1e-6; 

n = size(x,1); 
y = f(x); 
m = size(y,1); 
dy_dx = zeros(m, n); 

for i = 1:size(x,1)
    z = zeros(n,1); 
    z(i) = de; 
    x2 = x + z; 
    y2 = f(x2); 
    x1 = x - z; 
    y1 = f(x1);  
    dy_dx(:,i) = (y2-y1)/(2*de); 
    
end

end

