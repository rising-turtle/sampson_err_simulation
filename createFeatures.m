%%
% create a number of features by distributing around area [x, y, z] 
%   [-x/2, x/2], [-y/2, y/2] [0, z]
% 


function vfeats = createFeatures(x, y, z)

    vfeats = []; 
    id = 1;
    for i=-x/2:1:x/2
        for j = -y/2:1:y/2
            for k = 1:1:z
                feature.x = i;
                feature.y = j; 
                feature.z = k;
                feature.id = id; 
                id = id + 1;
                vfeats = [vfeats; feature];
            end
        end
    end

end