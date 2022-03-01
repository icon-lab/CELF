function rotData = rotateDataPoints(data,rotAngle)
%rotateDataPoints rotate points about origin by rotAngle(deg)
%   ----inputs----
%       data: (Nx2) data points
%       rotAngle: angle of rotation (clockwise) (radians)
%   ----output----
%       rotData: (Nx2) rotated data points

a = rotAngle;

% rotation matrix
R = [cos(a),-sin(a); sin(a),cos(a)];

% rotation about origin
rotData = (R*data')';

end