function [ u ] = windprofil( z, z0, m, zi, u0 )

%   windprofil calculates the wind speed at height z
%   Verbindung des logarithmischen Windprofils mit dem Potenzansatz, 
% da das logarithmische Windprofil nur in Prandtl-Schicht (~1/10 von zi) gueltig

% Eingabedaten
% z: Hoehe der Windgeschwindigkeit
% z0: Rauhigkeitslaenge
% zi: Hoehe der Grenzschicht
% u0: Referenzgeschwindigkeit in 10 m Hoehe
if z <= zi/10
    u = u0 * log(z/z0)/log(10/z0);
else
    ulog = u0 * log((zi/10)/z0)/log(10/z0);
    u = ulog * (z/(zi/10))^m;
end