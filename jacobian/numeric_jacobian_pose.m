%%
%    April 30, 2020, He Zhang, hzhang8@vcu.edu
%       numeric jacobians where y = f(x) is a matrix, 
%           x is a pose x = [t q]'
%           output dy_dx
%

function [dy_dx] = numeric_jacobian_pose(f,x)

if nargin <= 0
    % f = @(x)e2R(x); 
    % x = [1., 2.,0.5]';
    % [R, dR] = e2R(x);
end

de = 1e-6; 

% n = size(x,1);
n = 6; 
y = f(x); 
m = size(y,1); 
dy_dx = zeros(m, n); 

for i = 1:n
    z = zeros(n,1); 
    z(i) = de;
    x1 = x; 
    x2 = x; 
    if i <= 3
        x2(1:3) = x(1:3) + z(1:3);
    else
        dq2 = delta_q(z(4:6)); 
        q2 = qProd(x2(4:7), dq2); 
        x2(4:7) = q2; 
    end
    y2 = f(x2);
   
    if i <= 3
        x1(1:3) = x(1:3) - z(1:3); 
    else
        dq1 = delta_q(-z(4:6)); 
        q1 = qProd(x1(4:7), dq1); 
        x1(4:7) = q1; 
    end
    y1 = f(x1); 
    dy_dx(:,i) = (y2-y1)/(2*de); 
    
end



end

