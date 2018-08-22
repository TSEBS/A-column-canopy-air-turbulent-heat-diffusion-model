function [Rn,G0,H,LE,ustar,evap_fr,EF,EF0,varargout]= EnergyBalance_mapRn(d0, z0m, z0h,LCT,hc, fc,LAI, ..., 
                                                                            albedo,emissivity,SWd,LWd, LST_K,...
                                                                            hpbl, Zref, Tref_K, Uref, Earef,qaref, Pref, Ps,Gc)
%---------------------------------------------------------------------
% syntax: outputs [ ustar, LE, LE0, G0, Rn, H_DL, H_WL, H_i, evap_fr, re_i]= EnergyBalance( d0, z0m, z0h, fc, ..., 
%         inputs:[ albedo, emissivity, LST_K, SWd, LWd,hpbl, Zref, Tref_K, Uref, Earef, qaref, Pref, Ps, ri_i)
%---------------------------------------------------------------------
% Description of parameters
% Zref      -> Reference height                         (m)
% hpbl      -> Height of the PBL                        (m)
%              (if not available, use 1000m)
% d0        -> Zero plane displacement height           (m)
% z0m       -> Roughness height for momentum tranfer    (m)
% z0h       -> Roughness height for heat tranfer        (m)
% fc        -> Fractional vegetaion cover               (-)
% Uref      -> Wind speed at reference height           (m/s)
% Tref_K    -> Air temperature at reference height      (K)
% Pref      -> Pressure at reference height             (Pa)
% qaref     -> Specific humidity at reference height    (kg/kg)
% LST_K     -> Surface temperature                      (K)
% Ps        -> Surafce pressure                         (Pa)
% SWd       -> Downward Solar Radiation                 (Watt/m^2)
% LWd       -> Downward long wave radiation             (Watt/m^2)",
% albedo    -> Albedo                                   (-)
% emissivity-> Emissivity of the surface                (-)
%---------------------------------------------------------------------
% Solving of 3 equations
% Nots: Here we start to solve the system of three equations
% i.e. the equation of momentum transfer, the equation of heat transfer, and the equation for the stability length.
% We use ASL functions of Brutsaert, 1999, if Zref < hst, the ASL height
% We use BAS functions of Brutsaert, 1999, if Zref > hst, Zref <= hpbl, the
% PBL height.
% Note: We will determine the Monin-Obukhov length using the definition
% given by Brutsaert (1982), P.65, eqn. 5.25.
% i.e. L = - ustar^3*rhoa/(k*g*((H/Ta*Cp)+0.61*E))
% By using the quantity, ef=Le*E/(Rn-G0), we can write
% H = (1-ef)*(Rn-G0) and E = ef/Le*(Rn-G0)
% So that L =-ustar^3*rhoam/(k*g*((1-ef)/(Ta*CP)+0.61*ef/Le)*(Rn-G0))
% From this eqn., it is obvious that L=f(ustar^3) and L=f(ef^-1)
% LIMITING CASES: A good idea could be to take L_DL and L_WL respectively as
% the max and min stability length.
% This has the advantage that the feedback effects between land
% surface and the hpbl are considered.
% For theoretical limits, H=Rn-G0, and E=Rn-G0 respectively.
% (Note: For wet surfaces in a dry climate, E may be bigger than
% Rn-G0 due to negative H,
% i.e. the Oasis effect. In this case, a small modification of
% L_WL can be expected. This is ignored for the time being)
% Previously, we interpolated mu_i between mu_i(0) and
% mu_i(-1500), by means of temperature difference
% Though other factors may aslo have influences as seen in the
% calculation of T0Ta_l and T0Ta_u,
% the uncertainties associated to those factors do not warrant
% their adequate applications for this case.
% This is consistant with the definition of resistences for the
% limiting cases, that is, that we ignore the stable case
% over the whole region.
%% Constants
global P0 g k  Cpw Cpd L_e; 
global Rd Rv gamma Sigma_SB zs2hc phiuc Lambda_s Lambda_c; 
%--------------------added in 2017 July for parallization
P0      =   101325.;    % Standard pressure (Pa)
g       =   9.81;       % Gravity accelaration (kg s-2),[m/s^2]
k       =   0.4;        % von Karman constant
Cpw     =   1846;       % especific heat for water vapor, J Kg-1 K-1
Cpd     =   1005;       % especific heat for dry air, J Kg-1 K-1
L_e     =   2.430e+06 ; % J Kg-1
Rd      =   287.04;     % Gas Constant for Dry air, from table 2.1 P25 of Brutsaert 2005 (J kg-1 K-1)
Rv      =   461.5;      % Gas Constant for Water vapor, from table 2.1 P25 of Brutsaert 2005 (J kg-1 K-1)
gamma   =   67;         % psychometric constant (Pa K-1)
Sigma_SB=   5.678E-8;   % Stefan-Boltzmann's constant (W/m2/K4)
zs2hc   =   3.5;
phiuc   =   3;
%-----------------------------------Initialization
Rn                  =   LWd;
H                   =   LWd;
LE                  =   LWd;
G0                  =   LWd;
mlai                =   nanmax(LAI(:));
LAI(LCT==1)         =   mlai;  %
I                   =   LCT==1|LCT==2;
fc(I)               =   1-exp(-0.5*LAI(I));  % HETESSEL, page33,Massman Jour. Hyd 1999  
tau                 =   0;% 
beta                =   0.5;
fc                  =   1-exp(-beta.*(LAI./cos(tau))); % Eq. 4 in Murray and Verhoef 2007b
fc(LCT==4)          =   1-exp(-0.5*LAI(LCT==4));  % HETESSEL, page33,Massman Jour. Hyd 1999 
I                   =   LCT==1;
LAI(I)              =   LAI(I)./fc(I);  %
%-----------------------------------
Nu                  =   1.327e-5 * (P0./Pref).*((Tref_K ./ 273.15).^1.81);  % Kinematic viscosity of air (Massman 1999b) (10 cm^2/s)
Ps                  =   Pref.*exp(hc./(18400*(1+Tref_K./273.15))); % ѹ�߻��㹫ʽ�� Z2-Z1=18400��1+t/273��log( P1/P2); 
Lv                  =   (2.501-0.00234*(Tref_K-273.15))*1000000;   % Eq. 2.20 in Nicolas
%% Meteorological Parameters
Earef                               = 	Pref .* qaref * (Rv / Rd);         % actual vapour pressure (based on Pressure at reference height)
Theta_a                             =   Tref_K .*((P0./Pref ).^0.286);     % potential Air temperature    Eq. 2.23 P32 Brutsaert 2008
Theta_s                             =   LST_K  .*((P0./Ps   ).^0.286);     % potential surface temperature
Theta_av                            =   Theta_a.* (1 + 0.61 * qaref);      % virtual potential air temperature  P32 Brutsaert 2008 
% air densities (Eq 2.4,2.6 P25 Brutsaert 2008)
rhoa_m                              =   (Pref - 0.378 * Earef)./ (Rd * Tref_K);              % density of moist air [kg/m3]
rhoa_WL                             =   (Pref - 1.000 * Earef)./ (Rd * LST_K);               % density of dry air. [kg/m3]
% NOTE: rhoa_WL is only used for the wet-limit. To get a true upperlimit for the sensible heat
% the Landsurface Temperature is used as a proxy instead of air temperature.
% moist air density (kg / m3), this the same as rhoa_m
Cp                                  =   qaref* Cpw + (1-qaref)*Cpd;                          % Specific heat for constant pressure Page 29 of Brutsaert 2005
rhoa_m_Cp                           =   rhoa_m .* Cp;                                        % specific air heat capacity (J K-1 m-3)
rhoa_WL_Cp                          =   rhoa_WL .* Cp;                                       % specific air heat capacity (J K-1 m-3)
mw                                  =   0.622;     % ratio molecular weight of water vapor/dry air
gamma                               =   Cp.*Pref./(Lv*mw);    % psychrometric constant (Pa K^-1) when Pref in unit Pa 
% Computing saturated vapor pressure for land surface at wet limit (WL) (this is for lower limit of PBL)
LST_C                               =   LST_K - 273.15;                                     % Landsurface temperature[C].
A                                   =   611;                                                % [Pa]
B                                   =   17.502;
C                                   =   240.97;                                             % [C]
esat_WL                             =   A * exp(B * LST_C./(LST_C + C));                    % Pa,(3.8),p.41 of Campbell & Norman, 1998
slope_WL                            =   B * C * esat_WL ./ ((C + LST_C).^2);                % Pa*0C-1,(3.9)
% NOTE: esat_WL is only used for the wet-limit. To get a true upperlimit for the sensible heat
% the Landsurface Temperature is used as a proxy instead of air temperature.
%% Rn and G0
% Rn                                =   (1.0 - albedo) .* SWd + emissivity .* LWd - emissivity .* Sigma_SB .* LST_K.^4;
Rn                                  =   (1.0 - albedo) .* SWd +  LWd - emissivity .* Sigma_SB .* LST_K.^4;
Lambda_s                            =	0.315;     % bare soil (Kustas et al., 1989)
Lambda_c                            =	0.050;     % full vegetation canopy (Monteith, 1973)
G0                                  =   Rn .* (Lambda_c + (1 - fc) .* (Lambda_s - Lambda_c));
%%
G0(LCT==0)  = Rn(LCT==0).*0.5;      % water
G0(LCT==11) = Rn(LCT==11).*0.4;     % wetland
% %% snow and ice
G0(LCT==15) = Rn(LCT==15).*0.05; 
% %% Urban
G0(LCT==13) = Rn(LCT==13).*0.15; 
%% ASL height
% hst= alfa*hpbl, with alfa=0.10~0.15 over moderately rough surfaces, or hst=beta*z0, with beta= 100~150.
alfa                                =   0.12;                                               % These are mid values as given by Brutsaert,1999
beta                                =   125;                                                % These are mid values as given by Brutsaert,1999
hst                                 =   max(alfa * hpbl, beta * nanmax(z0m(:)));            % height of ASL
%% U* and L
CH                                  =   (Theta_s - Theta_a) .* k .* rhoa_m_Cp;              % from Eq. 2.55 P47 of Brutsaert 2008
CL                                  =   -rhoa_m_Cp .* Theta_av/ (k * g);                    % in this formula, air virtual potential temperature
z_d0                                =   Zref - d0;
ku                                  =   k * Uref;  
log_z_d0_z0m                        =   log(z_d0 ./ z0m);
log_z_d0_z0h                        =   log(z_d0 ./ z0h);
%% Forest sensible heat flux
if LCT(1,1)<=5&&LCT(1,1)>=1  % forest
    dt                           =   d0;  % 
    zstar                        =   zs2hc.*hc; 
    L                            =   -10^6;    % initial L 
    ustar                        =   Uref.*k./(log((Zref-dt)./z0m));    % equation 21 in Kitsiriw 2012
    H                            =   CH .* ustar ./ log((Zref-dt)./z0h);% H  in neutral condition when stability factors are zero
    %-------------------------H initial value
    I    =  (H-Rn)>10;
    H(I) =  0.52*Rn(I);   
    I    =  H<-20;
    H(I) =  0.52*Rn(I);    
    H(imag(H)~=0)= NaN;   
     %-------------------------
    H0                      =   H;                                              % H0 is H in neutral condition
    errorH                  =   10;
    steps                   =   0;
    while   (nanmax(errorH(:)) > 0.1 && steps < 100)  
        L                   =   CL .* (ustar.^3) ./ H;
        psu1                =   PSIu1(hc,LAI,Zref,zstar,d0,L,1);
        psim1               =   PSIm((Zref-dt)./L);
        psim2               =   PSIm(z0m./L);
        ustar               =   ku ./ ( log((Zref-dt)./z0m) - psim1 + psim2 + psu1); % equation 2 in Kitsri, 2013 AFM
        psu2                =   PSIu1(hc,LAI,Zref,zstar,d0,L,2);
        psih1               =   PSIh((Zref-dt)./L);
        psih2               =   PSIh(z0h./L);
        H                   =   CH .* ustar ./ ( log((Zref-dt)./z0h) - psih1 + psih2 + psu2);
        errorH              =   abs(H0 - H);
        H0                  =   H;
        steps               =   steps + 1;
    end
    C_i1                    =   PSIh(Zref./ L);
    C_i2                    =   PSIh(z0h ./ L);
else
% Initial guess for u*, H and L assuming neutral stability
L                           =   0;                                                  % initial L is zero for neutral condition
ustar                       =   ku ./ log_z_d0_z0m;                                 % U* in neutral condition when stability factors are zero
H                           =   CH .* ustar ./ log_z_d0_z0h;                        % H  in neutral condition when stability factors are zero
%-------------------------H initial value
I    =  (H-Rn)>10;
H(I) =  0.52*Rn(I);   
I    =  H<-20;
H(I) =  0.52*Rn(I);    
H(imag(H)~=0)= NaN;   
%-------------------------
errorH                              =   10;
H0                                  =   H;                                                  % H0 is H in neutral condition
steps                               =   0;
 if max(Zref(:)) <= hst
    while   (max(errorH(:)) > 0.1 && steps < 100) % MOS  
        L                           =   CL .* (ustar.^3) ./ H;
        psim1                       =   PSIm(z_d0./L);
        psim2                       =   PSIm(z0m./L);
        ustar                       =   ku ./ (log_z_d0_z0m - psim1 + psim2);
        psih1                       =   PSIh(z_d0./L);
        psih2                       =   PSIh(z0h./L);
        H                           =   CH .* ustar ./ (log_z_d0_z0h - psih1 + psih2);
        H(H>600)                    =   NaN;% 
        H(H<-50)                    =   NaN;
        errorH                      =   abs(H0 - H);
        H0                          =   H;
        steps                       =   steps + 1;        
    end
    C_i1                            =   PSIh(Zref./ L);
    C_i2                            =   PSIh(z0h ./ L);
else
   while    (max(errorH(:)) > 0.1 && steps < 100)            % BAS
        L                           =   CL .* (ustar.^3) ./ H;                              % this is temporary L computed based on previous step Ustar and H
        ustar                       =   ku ./ (log_z_d0_z0m - Bw(hpbl, L, z0m, d0));        % Eq. 2.67 P 52 Brutsaert 2008
        H                           =   CH .* ustar./(log_z_d0_z0h - Cw(hpbl,L,z0m,z0h,d0));% Eq. 2.68 P 52 Brutsaert 2008
        errorH                      =   abs(H0 - H);
        H0                          =   H;
        steps                       =   steps + 1;
   end
   C_i1                             =   Cw(Zref, L, z0m, z0h, d0);
   C_i2                             =   0;
 end
end
%-------------------------------------
I    =  (H-Rn)>10;
H(I) =  0.52*Rn(I);   %
I    =  H<-20;
H(I) =  0.52*Rn(I);    %
H(imag(H)~=0)= NaN;    %      
H(SWd<50&LCT==1)=NaN; % quality control for night time,new H parameterization doesnt work on night
H(isnan(SWd))=NaN;    % quality control for night time,new H parameterization doesnt work on night
%% Sensible heat Flux
I_1                                 =   log_z_d0_z0h + C_i2 >  C_i1;
I_2                                 =   log_z_d0_z0h + C_i2 <= C_i1;

re_i1                               =   (log_z_d0_z0h - C_i1 + C_i2)./(k * ustar);          % Actual resistance to heat transfer (s/m)
re_i2                               =   (log_z_d0_z0h      )./(k * ustar);
re_i                                =   I_1.*re_i1 + I_2.*re_i2;

H_i                                 =   rhoa_m_Cp.* (Theta_s-Theta_a)./re_i;                % Sensible heat flux
%% Sensible heat flux at theoretical Dry limit
% Dry limit
% L_dry                              = 	0;
H_DL                               =   Rn - G0;
%% Sensible heat flux at theoretical wet limit
% Dry air is assumed.. eact=0, and we need to take the density of dry air
L_e                                 =   2.501-0.002361*(Tref_K-273.15); % latent heat of vaporation (MJ kg-1), 
L_e                                 =   L_e*1000000; % transfer to :J kg-1
L_WL                                =   -(ustar.^3).* rhoa_WL./(k*g*(0.61* (Rn - G0)./ L_e)); 
% Bulk Stability Corrections
I_MOS                               =   (Zref < hst);
I_BAS                               =   (Zref >= hst);
C_WL_MOS                            =   PSIh(Zref ./ L_WL); %  minus should be plus, according to Joris.
C_WL_BAS                            =   Cw(Zref, L_WL, z0m, z0h,d0);
C_WL                                =   I_MOS.*C_WL_MOS + I_BAS.*C_WL_BAS;
% Calculating Resistances
I_1                                 =   log_z_d0_z0h>C_WL;
I_2                                 =   log_z_d0_z0h<=C_WL;
re_WL1                              =   (log_z_d0_z0h - C_WL)./(k * ustar);                 % Actual resistance to heat transfer at wet limit(s/m)
re_WL2                              =   (log_z_d0_z0h       )./(k * ustar);
re_WL                               =   I_1.*re_WL1 + I_2.*re_WL2;
H_WL                                =   ((Rn - G0) - (rhoa_WL_Cp./re_WL).*((esat_WL)./gamma))./(1.0 + slope_WL./gamma);
%% Evaporative fraction 
% H_i                                 =   min(H_i, H_DL);
H_i(H_i>H_DL)                       =   H_DL(H_i>H_DL); % 
% H_i                                 =   max(H_i, H_WL);
H_i(H_i<H_WL)                       =   H_WL(H_i<H_WL); % 
% Relative evaporation
I_water                             =   (H_DL <= H_WL);
I_land                              =   (H_DL > H_WL);
evap_re_water                       =   1;                                                  % for water & wet surfaces
evap_re_land                        =   1 - (H_i - H_WL) ./ (H_DL - H_WL);
evap_re                             =   I_water.*evap_re_water + I_land.*evap_re_land;
% Evaporative fraction
I_1                                 =   ((Rn - G0) ~= 0);
I_2                                 =   ((Rn - G0) == 0);
evap_fr_1                           =   evap_re .* (Rn - G0 - H_WL) ./ (Rn - G0);
evap_fr_2                           =   1;                                                  % for negative available energy
evap_fr                             =   min(I_1.*evap_fr_1 + I_2.*evap_fr_2,1);
%% Latent heat
LE                                  =   evap_fr .* (Rn - G0);                                % Latent heat flux W m-2
LE(imag(LE)~=0)                     =   NaN;    % 

%LE                                  =   Rn-G0-H;
return;
