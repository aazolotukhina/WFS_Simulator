function [WF_derivX,WF_derivY] = zernike_8_deriv(X,Y,z11,zi11,z02,z22,zi22,z13,zi13,z04)
% forms a wavefront derivatives by zernike polynomials

% Partial derivatives by X of Zernike polynomials defined for Cartesian 
% coordinates (x, y)
% INPUT:
% X,Y - coords of sensor
% z11 - Tip (X-Tilt, horizontal tilt)
% zi11 - Tilt (Y-Tilt, vertical tilt)
% z02 - Defocus (longitudinal position)
% z22 - Vertical astigmatism (0 deg)
% zi22 - Oblique astigmatism (45 deg)
% z13 - Horizontal coma
% zi13 - Vertical coma
% z04 - Primary spherical

% OUTPUT:
% WF_deriv - wavefront derivatives


% Zernike polynomials
WF_derivX = z11*1 + zi11*0 +...
    z02*4*X + z22*2*X + zi22*2*Y +...
    z13*(3*(3*X.^2+Y.^2)-2) + zi13*6*X.*Y +...
    z04*12*X.*(2*(X.^2+Y.^2)-1);

WF_derivY = z11*0 + zi11*1 +...
    z02*4*Y + z22*(-2)*Y + zi22*2*X +...
    z13*6*X.*Y + zi13*(3*(X.^2+3*Y.^2)-2) +...
    z04*12*Y.*(2*(X.^2+Y.^2)-1);
end

