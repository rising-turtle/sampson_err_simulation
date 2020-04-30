
function M = skew_symmetric33(x)
    
    M = [0 -x(3) x(2); 
         x(3) 0 -x(1);
         -x(2) x(1) 0];

end