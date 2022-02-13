function [WF_derivX,WF_derivY] = deriv(X,Y,R)
%forms a wavefront by sphere
%Partial derivatives by X and by Y
%INPUT
% X,Y - coords of sensor
% R - distance between the point source and the MLA-plane

% Derived spheres
WF_derivX = X./sqrt(R^2-X.^2-Y.^2);
WF_derivY = Y./sqrt(R^2-X.^2-Y.^2);

end

