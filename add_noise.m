
function obs_n = add_noise(obs, sigma)

    obs_n = obs; 
    for i = 1:size(obs,1)
        obs_n(i).pi_n.x = obs_n(i).pi_n.x + randn(1,1)*sigma; 
        obs_n(i).pi_n.y = obs_n(i).pi_n.y + randn(1,1)*sigma; 
        obs_n(i).pj_n.x = obs_n(i).pj_n.x + randn(1,1)*sigma; 
        obs_n(i).pj_n.y = obs_n(i).pj_n.y + randn(1,1)*sigma; 
    end

end
