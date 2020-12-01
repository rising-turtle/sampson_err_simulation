%%
% Nov. 30 2020, He Zhang, fuyinzh@gmail.com
% 
% triangulate point given measurement pt0, pt1 observed at Pose0, Pose1
%

function pt3d = triangulate_point(Pose0, Pose1, pt0, pt1)

    M = zeros(4); 
    M(1,:) = pt0(1)*Pose0(3, :) - Pose0(1, :); 
    M(2,:) = pt0(2)*Pose0(3, :) - Pose0(2, :); 
    M(3,:) = pt1(1)*Pose1(3, :) - Pose1(1, :); 
    M(4,:) = pt1(2)*Pose1(3, :) - Pose1(2, :); 
    
    [U, D, V] = svd(M); 
    v4 = V(:,4); 
    pt3d = zeros(3,1);
    pt3d(1) = v4(1)/v4(4); 
    pt3d(2) = v4(2)/v4(4); 
    pt3d(3) = v4(3)/v4(4); 

end