function [ u ] = windprofil( z, z0, m, zi, u0 )

%   windprofil calculates the wind speed at height z
%   Verbindung des logarithmischen Windprofils mit dem Potenzansatz, 
% da das logarithmische Windprofil nur in Prandtl-Schicht (~1/10 von zi) gueltig

% Input values
% z: Height of windspeed
% z0: Roughness length
% zi: Height of boundary layer
% u0: wind speed in height of 10 m
if z <= zi/10
    u = u0 * log(z/z0)/log(10/z0);
else
    ulog = u0 * log((zi/10)/z0)/log(10/z0);
    u = ulog * (z/(zi/10))^m;
end
