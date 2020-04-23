function obs = overlap_obs(obs_i, obs_j)
Ni = size(obs_i, 1); 
Nj = size(obs_j, 1);
obs = []; 
for i =1:Ni
    id_i = obs_i(i).feat_id;
    good = 0; 
    jj = -1;
    for j = 1:Nj
        id_j = obs_j(j).feat_id; 
        if id_i == id_j
            jj = j; 
            good = 1; 
            break; 
        end
    end
    
    if good == 1 % find one 
        obs_ij.pi.x = obs_i(i).obs_x; 
        obs_ij.pi.y = obs_i(i).obs_y; 
        obs_ij.di = obs_i(i).lpt(3); %% depth 
        obs_ij.pj.x = obs_j(jj).obs_x; 
        obs_ij.pj.y = obs_j(jj).obs_y; 
        obs_ij.pi_n = obs_ij.pi; 
        obs_ij.pj_n = obs_ij.pj;
        obs = [obs_ij; obs];
    end
    
end