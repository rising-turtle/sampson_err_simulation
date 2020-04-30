%% update quaternion q = [q, b, c, d] = [qw, qx, qy, qz]
function q = delta_q(theta)
    
    q = zeros(4,1); 
    half_theta = theta./2; 
    q(1) = 1.; 
    q(2) = half_theta(1); 
    q(3) = half_theta(2); 
    q(4) = half_theta(3); 

end