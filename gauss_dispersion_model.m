%% Gauss dispersion model

% This model calculates the dispersion of a polution cloud coming out of a chimney
% I created this model for a project with the the departement of technical meteorology at the University of Hamburg

clear; close all; clc
%% Declaration
% fixed parameters
zi = 800;   % height of boundary layer [m]
x  = 0:1:2000 % dispersial area up to 2km [m]
a  = 1; % Deposition (is neglected)

% Output parameters
H_s = 65;   % heigt of chimney
u0_s = 4.5; % wind speed [m/s]
m_s = 0.31; % wind profile exponent for neutral Stratification
            % 0.31 for urban roughness; 0,18 for smooth roughness

% roughness
z0_s = 2; % roughness length [m]
Q_s = 1; % stream of emissions [kg/h]

% varied parameters
H = [40 50 H_s 80 90]; % height of chimney [m]
u0 = [1.5 3 u0_s 6 7.5 9]; % wind speed [m/s]
m = [0.18 m_s]; % wind profile exponent
z0 = [0.3 z0_s]; % roughness length
Q = [Q_s 25 50 75 100]; % stream of emissions

% Function to calculate the wind profile
u_s = windprofil(H_s, z0_s, m_s, zi, u0_s);

%% wind profile

z = 10:1:100;

u_1 = zeros(1,length(z));
u_2 = zeros(1,length(z));

for i = 1:length(z)
    u_1(i) = windprofil(z(i), z0(1), m(1), zi, u0_s);
    u_2(i) = windprofil(z(i), z0(2), m(2), zi, u0_s);
end

figure(5);
plot(u_1, z, 'linewidth', 1.5);
hold on;
plot(u_2, z, 'linewidth', 1.5);
hold on;
plot([0 12], [zi/10 zi/10], 'linestyle', '--', 'col', 'k');
grid on;
legend('z0 = 0.3m', 'z0 = 2m', 'location', 'northwest');
xlabel('wind speed u (ms^{-1})');
ylabel('height z (m)');


%% Calculate SIGMA values

% parameters for H < 50m
F50 = 0.640;
f50 = 0.784;
G50 = 0.215;
g50 = 0.885;

% parameters for H = 100m
F100 = 0.504;
f100 = 0.818;
G100 = 0.265;
g100 = 0.818;

% linear and logarithmic fitting of dispersion parameters for dispersion height in between
if H_s < 50
    F = F50;
    f = f50;
    G = G50;
    g = g50;
elseif H_s >= 50 && H_s < 100
    % log. fit
    F = F50 * exp((H_s - 50)* (log(F100) - log(F50)) / (100-50));
    G = G50 * exp((H_s - 50)* (log(G100) - log(G50)) / (100-50));
    % lin. fit
    f = f50 + (f100 - f50) / (100 - 50) * (H_s - 50);
    g = g50 + (g100 - g50) / (100 - 50) * (H_s - 50);
end

sigma_y = F * x.^f;
sigma_z = G * x.^g;


%% calculating concentration
% with variation of dispersion parameters

% 1. Variation of dispersion height (H)

C0_variH = zeros(5,length(x));

for i = 1:length(H)
    for j = 1:length(x)
        u = windprofil(H(i), z0_s, m_s, zi, u0_s);
        C0_variH(i,j) = ((Q_s/3600) / (u*2*pi*sigma_y(j)*sigma_z(j))) * (exp(-0.5*((0-H(i))/sigma_z(j)).^2)+a*exp(-0.5*((0+H(i))/sigma_z(j)).^2));
    end
end

figure(1);
plot(x, C0_variH, 'linewidth', 1.5);
legend('H = 40', 'H = 50', 'H = 65', 'H = 80', 'H = 90');
grid on;
ylim([0 3*10^(-9)]);
xlabel('Distance x (m)');
ylabel('normalized concentration C (1)');


% 2. Variation of wind speed

C0_variu0 = zeros(6,length(x));

for i = 1:length(u0)
    for j = 1:length(X)
        u = windprofil(H_s, z0_s, m_s, zi, u0(i));
        C0_variR(i,j) = ((Q_s/3600) / (u*2*pi*sigma_y(j)*sigma_z(j))) * (exp(-0.5*((0-H_s)/sigma_z(j)).^2) + a*exp(-0.5*((0+H_s)/sigma_z(j)).^2));
    end
end

figure(3);
plot(x,C0_variR,'linewidth',1.5);
legend('städtische Umgebung','ländliche Umgebung');
grid on;
ylim([0 3*10^(-9)]);
xlabel('Entfernung x (m)');
ylabel('normierte Konzentration C (1)');


% 3. Variation der Rauigkeit

C0_variR = zeros(2,length(x));

for i = 1:length(z0)
    for j = 1:length(x)
        u = windprofil( H_s, z0(i), m(i), zi, u0_s );
        C0_variR(i,j) = ((Q_s/3600)/(u*2*pi*sigma_y(j)*sigma_z(j)))*(exp(-0.5*((0-H_s)/sigma_z(j)).^2)+a*exp(-0.5*((0+H_s)/sigma_z(j)).^2));
    end
end

figure(3);
plot(x,C0_variR,'linewidth',1.5);
legend('urban environment','rural environment');
grid on;
ylim([0 3*10^(-9)]);
xlabel('Distance x (m)');
ylabel('normalized Concentration C (1)');


% 4. Variation of stream of emission

C0_variQ = zeros(5,length(x));

for i = 1:length(Q)
     for j = 1:length(x)
        C0_variQ(i,j) = ((Q(i)*10^6/3600)/(u_s*2*pi*sigma_y(j)*sigma_z(j)))*(exp(-0.5*((0-H_s)/sigma_z(j)).^2)+a*exp(-0.5*((0+H_s)/sigma_z(j)).^2));
     end
end

figure(4);
 plot(x,C0_variQ,'linewidth',1.5);
 legend('Q = 1 kg\cdot{}h^{-1}','Q = 25 kg\cdot{}h^{-1}','Q = 50
 kg\cdot{}h^{-1}','Q = 75 kg\cdot{}h^{-1}','Q = 100 kg\cdot{}
 h^{-1}');
 grid on;
 xlabel('Distance x (m)');
 ylabel('Concentration C (mg\cdot{}m^{-3})');


%% Speichern der Plots
%print(figure(1),'/home/zmaw/u231124/Technische Meteorologie/
%GaussBodenkonzentration_VariationSchornsteinhoehe','-dpng');
%print(figure(2),'/home/zmaw/u231124/Technische Meteorologie/
%GaussBodenkonzentration_VariationWindgeschwindigkeit','-dpng');
%print(figure(3),'/home/zmaw/u231124/Technische Meteorologie/
%GaussBodenkonzentration_VariationRauhigkeit','-dpng');
%print(figure(4),'/home/zmaw/u231124/Technische Meteorologie/
%GaussBodenkonzentration_VariationEmissin','-dpng');
%print(figure(5),'/home/zmaw/u231124/Technische Meteorologie/
%Windprofil_StadtLand','-dpng');
   