%%
% test whether the linearity is still behold 
% 
% It seems like that the sampson error itself is a linear function, 
% linearity test always get %5 pass rate, no meaning to do this
%
%

function rej_array = linearity_check_for_sampson_err(N, sig_array)

    if nargin <= 0 
        N = 10; 
        sig_array = [1:0.2:2]';
    end
    

    rej_array = [] ;
    
    for i = 1:size(sig_array,1)
        
        rej_rate = monte_carlo_test(N, sig_array(i)); 
        
        rej_array = [rej_array; rej_rate]; 
        
        fprintf('rej_rate = %f noise_sigma = %f \n', rej_rate, sig_array(i));
        
    end
    

end

function rej_rate = monte_carlo_test(N, noise_sigma)

%% generate features in a range
feats = createFeatures(8, 8, 7);

%% create observations for the first camera pose 
R = [1 0 0; 0 1 0; 0 0 1]; 
t = [0 0 0]'; 
cam = get_struct_core(); 
obs_i = createObservations(feats, cam, R, t); 

%% number of trials that pass chi2 test 
good = 0; 
bad = 0; 

%% random 
for k = 0:N
    
    % randomly generate next camera pose within rpy [-30, 30] degree, xyz
    rr = randn(6,1); 
    s = max(rr); 
    if s < -min(rr)
        s = -min(rr);
    end
    rr(1:3) = rr(1:3)/s*30; 
    rr(4:6) = rr(4:6)/s*2; 
    euler_angle = rr(1:3)*pi/180;
    R = e2R(euler_angle);
    t = rr(4:6); 
    %% find observations 
    obs_j = createObservations(feats, cam, R, t);
    
    %% overlap observations
    obs = overlap_obs(obs_i, obs_j); 
    if size(obs,1) < 10
        fprintf('too small observations for k = %f with %d obs', k+1, size(obs,1));
        continue; 
    end
    
    %% add noise
    obs = add_noise(obs, noise_sigma); 
    
    for i=1:10
       ne = compute_sampson_err(obs(i), cam, R, t, noise_sigma); 
       
       pass = chi2_test(ne); 
       
       if pass == 1
           good = good + 1; 
       else
           bad = bad + 1;
       end
    end
end

rej_rate = bad/(good+bad); 

fprintf('given noise = %f, reject rate is = %f with %d tests \r\n', noise_sigma, rej_rate, (good+bad));

end

function good = chi2_test(x)

    chi2 = x(1)*x(1) + x(2)*x(2);
    good = 0; 
    if chi2 <= 5.99
        good = 1;
    end

end

function ne = compute_sampson_err(obs, cam, R, t, sig)

    %% compute error
    nxi = (obs.pi_n.x-cam.cx)/cam.fx; 
    nyi = (obs.pi_n.y-cam.cy)/cam.fy; 
    nxj = (obs.pj_n.x-cam.cx)/cam.fx; 
    nyj = (obs.pj_n.y-cam.cy)/cam.fy; 
    d = obs.di; 
    x = [nxi nyi nxj nyj]'; 
    e = [0 0]';
    zj = (R(3,1)*x(1)*d + R(3,2)*x(2)*d + R(3,3)*d + t(3));
    e(1) = -x(4)*zj+(R(2,1)*x(1)*d + R(2,2)*x(2)*d + R(2,3)*d + t(2));
    e(2) = x(3)*zj-(R(1,1)*x(1)*d + R(1,2)*x(2)*d + R(1,3)*d + t(1)); 

    %% compute covariance matrix 
    de_dpi = [0 1 -nyj; -1 0 nxj]*R*[d 0; 0 d; 0 0]; 
    de_dpj = [0 -zj; zj 0]; 
    J = [de_dpi de_dpj];
    
    COV_X = eye(4)*(sig/cam.fx)^2; 
    COV_e = J*COV_X*J';
    [U,D, ] = svd(COV_e); 
    B = U*sqrt(D); 
    
    %% normalized 
    ne = inv(B)*e; 
    
end