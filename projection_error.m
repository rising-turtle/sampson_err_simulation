
function c = projection_error(R, t, x, d)
    
    pti = [x(1)*d, x(2)*d, d]'; 
    ptj = R * pti + t; 
    
    c1 = (ptj(3)*x(3) - ptj(1)); 
    c2 = (ptj(3)*x(4) - ptj(2)); 

    c = sqrt(c1^2+c2^2);
    
end