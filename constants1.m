%% Constants
global P0 T0 Rd Rv g k Sigma_SB Cpw Cpd gamma Pr L_e 
T0      =   273.15;     % zero Kelvin [C]
P0      =   101325.;    % Standard pressure (Pa)
Rd      =   287.04;     % Gas Constant for Dry air, from table 2.1 P25 of Brutsaert 2005 (J kg-1 K-1)
Rv      =   461.5;      % Gas Constant for Water vapor, from table 2.1 P25 of Brutsaert 2005 (J kg-1 K-1)
g       =   9.81;       % Gravity accelaration (kg s-2)
k       =   0.4;        % von Karman constant
Sigma_SB=   5.678E-8;   % Stefan-Boltzmann's constant (W/m2/K4)
Cpw     =   1846;       % especific heat for water vapor, J Kg-1 K-1
Cpd     =   1005;       % especific heat for dry air, J Kg-1 K-1

gamma   =   67;         % psychometric constant (Pa K-1)
Pr      =   0.7;        % Prandtl number

% ri_i    =   60;         % surface resistance of standard crop, s m-1,

% The latent heat of vaporization at 30C from Brutsaert, 1982, p.41,tab. 3.4,
% more exact values can be obtained from eqn(3.22, 3.24a,b)
L_e     =   2.430e+06 ; % J Kg-1