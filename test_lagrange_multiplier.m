% syms l b h lambda
% A = l * b + 2 * b * h + 2 * l * h; % area of needed sheet metal (thickness = 1, because we are lazy guys)
% V = l * b * h - 10 == 0; % volume constraint 10 m³
% L = A + lambda * lhs(V); % Langrange function
% dL_dl = diff(L,l) == 0; % derivative of L with respect to l
% dL_db = diff(L,b) == 0; % derivative of L with respect to b
% dL_dh = diff(L,h) == 0; % derivative of L with respect to h
% dL_dlambda = diff(L,lambda) == 0; % derivative of L with respect to lambda
% system = [dL_dl; dL_db; dL_dh; dL_dlambda]; % build the system of equations
% [l_val, b_val, h_val,lambda_val] = solve(system, [l b h lambda], 'Real', true) % solve the system of equations and display the results 
% results_numeric = double([l_val, b_val, h_val,lambda_val]) % show results in a vector of data type double
% V_res = l_val * b_val * h_val % calculate the resulting volume (should be 10, if correct)
% A_res = l_val * b_val + 2 * b_val * h_val + 2 * l_val * h_val % calculate the needed amount (area) of sheet metal

syms x y lambda 
F = (x-1.5)^2 + (y-0.3)^2; 
G = x^2 + y^2 - 1 == 0; 
L = F + lambda * lhs(G); 
dL_dx = diff(L, x) == 0; 
dL_dy = diff(L, y) == 0; 
dL_dlambda = diff(L, lambda) == 0; 
system = [dL_dx; dL_dy; dL_dlambda]; 
[lx, ly, la] = solve(system, [x y lambda], 'Real', true); 
result_numeric = double([lx, ly, la]); 

VF = zeros(size(lx, 1), 1);
VG = zeros(size(lx, 1), 1);
for i=1:size(lx)
    VG(i) = lx(i)^2 + ly(i)^2;
    VF(i) = (lx(i)-1.5)^2 + (ly(i)-0.3)^2;
end
VG
VF 
min(VF)


