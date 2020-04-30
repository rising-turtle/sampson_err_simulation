function [q,Qe] = e2q(e)

% E2Q Euler angles to quaternion conversion.
%   Q = vE2Q(E) gives the quaternion Q corresponding to the Euler angles
%   vector E = [roll;pitch;yaw].
%
%   [Q,Q_e] = E2Q(E) returns also the Jacobian matrix Q_e = dQ/dE.
%
%   See also QUATERNION, EULERANGLES, R2Q, Q2E, Q2V.

%   Copyright 2008-2009 Joan Sola @ LAAS-CNRS.


qx = au2q(e(1),[1;0;0]); %  roll rot. on X axis
qy = au2q(e(2),[0;1;0]); % pitch rot. on Y axis
qz = au2q(e(3),[0;0;1]); %   yaw rot. on Z axis

q = qProd(qProd(qz,qy),qx);

if nargout == 2

    sr = sin(e(1)/2);
    sp = sin(e(2)/2);
    sy = sin(e(3)/2);

    cr = cos(e(1)/2);
    cp = cos(e(2)/2);
    cy = cos(e(3)/2);

    Qe = 0.5*[
        [ -cy*cp*sr+sy*sp*cr, -cy*sp*cr+sy*cp*sr, -sy*cp*cr+cy*sp*sr]
        [  cy*cp*cr+sy*sp*sr, -cy*sp*sr-sy*cp*cr, -sy*cp*sr-cy*sp*cr]
        [ -cy*sp*sr+sy*cp*cr,  cy*cp*cr-sy*sp*sr, -sy*sp*cr+cy*cp*sr]
        [ -sy*cp*sr-cy*sp*cr, -cy*cp*sr-sy*sp*cr,  cy*cp*cr+sy*sp*sr]
        ];
end

return

%%
syms r p y real
e = [r;p;y];

[q,Qe] = e2q(e);

simplify(Qe-jacobian(q,e))



% ========== End of function - Start GPL license ==========


%   # START GPL LICENSE

%---------------------------------------------------------------------
%
%   This file is part of SLAMTB, a SLAM toolbox for Matlab.
%
%   SLAMTB is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%
%   SLAMTB is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with SLAMTB.  If not, see <http://www.gnu.org/licenses/>.
%
%---------------------------------------------------------------------

%   SLAMTB is Copyright 2007,2008,2009 
%   by Joan Sola, David Marquez and Jean Marie Codol @ LAAS-CNRS.
%   See on top of this file for its particular copyright.

%   # END GPL LICENSE

