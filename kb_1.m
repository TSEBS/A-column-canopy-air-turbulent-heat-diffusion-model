function [z0m, d0, z0h,kB_1,kB1_s, kB1_c,kB1_m,zeth, varargout]=  kb_1(yy,LCT,fc,NDVI, LAI, hc, Zref, Uref, Pref, Tref_K,LST_K,qaref,kbc,kbs,varargin)
% KB_1M by Massman, 1999, (24) This a surrogate for the full LNF model for describing the combined
% canopy and soil boundary layer effects on kBH_1. The formulation follows Su (2001), HESS
% zref, reference height (2hc <= zref <= hst=(0.1 ~0.15)hi, (=za))
% Uref: wind speed at reference height
% u_h: wind speed at canopy height
% hc: canopy height
% LAI: canopy total leaf area index
% Cd: foliage drag coefficient, = 0.2
% Ct: heat transfer coefficient
% fc: Fractional canopy cover (-)
% Ta: ambient air temperature (0C)
% Pa: air pressure (Pa)
% hs: height of soil roughness obstacles
% Nu: Kinematic viscosity of air, (m^2/s)
% Ta: ambient temperature of the air, (0C)
% pa: ambient pressure of the air, (0C)
global c1;global c2;global c3;global Ct;global hs1;global Cd;global fh;global A1;global A2;global A3;global hcp;global stp;global lmdars;
%% Constants
% Cd      =   0.2;      % Foliage drag coefficient
% Ct      =   0.01;     % Heat transfer coefficient,Binbin discusssion 0.1-0.001;,default 0.01
Pr      =   0.71;       % Prandtl number
k       =   0.4;        % von Karman constant
T0      =   273.15;     % zero Kelvin [C]
P0      =   101325.;    % Standard pressure (Pa)
Lf      =   1.5;        % the size of the leaves
g       =   9.81;       % Gravity accelaration (kg s-2)
[m,n]   =   size(fc);
hs(1:m,1:n)  =   hs1;   % momentum roughness parameter (0.009 ~ 0.024)(Su et al., 1997, IJRS, p.2105-2124.),default value is 0.0012
Nu                    =   1.327e-5 * (P0 ./ Pref) .* ((Tref_K / T0).^1.81);  % Kinematic viscosity of air (Massman 1999b) (10 cm^2/s)
%% U*/U(h)
fs                    =   1 - fc;
%% KB-1 for canopy only
if kbc==1  % SEBS 
    ust2u_h             =   c1 - c2 .* exp(-c3 .* Cd .* LAI) ;  % Cd*LAI =zeta(h),h: canopy height, Ratio of ustar and u(h) (Su. 2001 Eq. 8)
    Cd_bs               =   2*ust2u_h.^2;                    % Bulk surface drag cofficient (Su 2001)
    %% within-canopy wind speed profile extinction coefficient
    n_ec                =   Cd * LAI ./ (Cd_bs);  % Cd*LAI =zeta(h),h: canopy height, n_ec=zeta(h)/(2.*(ustar2u_h(h))^2.)
    I_1                 =   (n_ec ~= 0);
    I_2                 =   (n_ec == 0);
    kB1_c1              =   k * Cd ./(4 * Ct .* ust2u_h .* (1 - exp(-n_ec/2)));    % KB-1 for Full canopy only (Choudhury and Monteith, 1988)
    kB1_c2              =   0;                                                     % KB-1 for Full canopy only lower limit
    kB1_c               =   I_1.*kB1_c1 + I_2.*kB1_c2;                             % KB-1 for Full canopy only
    d2h                 =   1 - 1./(2*n_ec) .* (1 - exp(-2 * n_ec));               % Ratio of displacement height and canopy height (derived from Su 2002, eq 9)
    z0m                 =   hc.*(1-d2h).*exp(-k*ust2u_h.^(-1));   % roughness height for momentum (Eq 10 Su. 2001)
    d0                  =   d2h.*hc;             % displacement height
elseif kbc==5  % Chen et al. 2018 JGR
    %__________________________________________________________________________
    % Define z/h vector from top to bottom of the canopy: [0 <= z01 <= 1]
    %__________________________________________________________________________
    z0h  = 0.0025;
    delh = 0.0001;
    z01  = 0:delh:1;
    nz   = length(z01);
    uz0  = log(z01/z0h);
    vkar = 0.4;
    c1   = 0.38;
    c3   = 15;
    c2   = c1+vkar/log(z0h);
    [hazn]    = NormfoldisLAI(LCT(1)); % eq 2 in Massman 2018
    intbetam  = trapz(hazn);
    haz       = LAI*hazn/intbetam; % 
    [zetaz,hacpz,hacpn]   = Massdragdis(haz);
    canopywind            = 'hypercos'; % 'exponent';
    [ucan,usuh,nzet,nexp] = MassUProfileLAI(hacpn,nz,zetaz,canopywind);
    [ubc,usb,I1,doh,z0oh,rough,lroug] = ...
        calUUSLAII(ucan,uz0,z01,z0h,delh,vkar,usuh,nz,nzet,nexp,canopywind);
    d0                      =  doh*hc;
    z0m                     =  z0oh*hc;
    %----------------------------------------------------
    ust2u_h                 =   usuh;
    Cd_bs                   =   2*ust2u_h.^2;                   % Bulk surface drag cofficient (Su 2001),ust2u_h=ustar(h)/u(h),friction velocity at canopy height.
    % within-canopy wind speed profile extinction coefficient
    zeth                    =   hacpn(nz);
    zetz2zeth               =   zetaz;
    n_ec                    =   zeth ./ (Cd_bs);                    % equation 6 in Su et al. 2001; Cd*LAI =zeta(h),h: canopy height, n_ec=zeta(h)/(2.*(ustar2u_h(h))^2.)
    %--------------------Brusaert 1979 
    CL             =   sqrt(ust2u_h); % on page 368 in Brusaert 1979, Boundary-Layer Meteorology 16(4): 365-388.
    u_h            =   max(zeros(size(z0m)),Uref.* log((hc-d0)./z0m)./log((Zref-d0)./z0m));  % (within canopy) Horizontal windspeed at the top of the canopy
    u_z            =   u_h.*exp(-n_ec.*(1-zetz2zeth));
    Ctz            =   nanmean(u_z);
    Ctz            =   CL.*Pr^(-0.67).*(Ctz.*hcp./Nu).^(-0.5);% on page 368 in Brusaert 1979, Boundary-Layer Meteorology 16(4): 365-388.
%     Ctz            =   0.01; % SEBS
    kB1_c          =   k * Cd ./(40 * Ctz .* ust2u_h .* (1 - exp(-n_ec/2)));    % KB-1 for Full canopy only (Choudhury and Monteith, 1988)
end
%% KB-1 for mixed canopy and soil
I_1                 =   Nu~= 0; % classify pixels (mixed canopy and soil)
u_h                 =   max(zeros(size(z0m)),Uref.* log((hc-d0)./z0m)./log((Zref-d0)./z0m));  % (within canopy) Horizontal windspeed at the top of the canopy
Ustar_m             =   ust2u_h.*u_h;                       % friction velocity
Re_m                =   zeros(size(I_1));                   % roughness Reynolds number lower limit, Su et al. 2001
Re_m(I_1)           =   hs(I_1).*Ustar_m(I_1)./Nu(I_1);     % hs being the roughness height of the soil,Su et al. 2001
Ct_m                =   Pr^(-2/3) * (Re_m.^(-1/2));         % heat transfer coefficient for soil, Su et al. 2001
% Ct_m          = Ct;  % 2015, new method to make kb_1m small for high canopy
kB1_m                 =   (k * ust2u_h).* (z0m./ hc)./ Ct_m;   % KB-1 for Mixed canopy (Choudhury and Monteith, 1988)
%% KB-1 for soil only
% Ustar_s             =   Uref * k ./ log(Zref ./ hs);    % Friction velocity in case of soil only. (Brutsaert 2008, Eq 2.41 P43 )[m/2]
% I_1                 =   (Nu ~= 0 );
% I_2                 =   (Nu == 0 );
% Re_s1               =   hs.* Ustar_s ./ Nu;             % roughness Reynolds number, = hu*/v,Su et al. 2001
% Re_s2               =   0.0;                            % roughness Reynolds number, = hu*/v,Su et al. 2001
% Re_s                =   I_1.*Re_s1 + I_2.*Re_s2;        % roughness Reynolds number, = hu*/v,Su et al. 2001
Re_s                  =   Re_m;                           % Su et al. 2002
if kbs==1
    kB1_s               =   2.46 * (Re_s.^(1/4)) - log(7.4);     % KB-1 for Soil only (Brutsaert,1982),BR82
elseif kbs==5
    %     Zref=Zref(1,1);  % When Zref is matrix, yang function accept a single Zref value
    z0m1                =   z0m;   %  intialize the matrix
    Zref1               =   10;    % 
    z0m1(:)             =   hs1;%  
    [z0h ]              =   yang_kb_1(z0m1,Zref1,Zref1,Uref,LST_K+hc*0.0098,Tref_K,qaref, Pref);  % Yang et al (2002)
    z0h(isnan(LST_K))   =   NaN; % 2016, April, discovered that when LST_K is a nan value,z0h is not nan;
    kB1_s               =   log(z0m1./z0h);
end
%% KB-1 Combined all
% oumega    =
% fc        =  1-exp(0.5*oumega*LAI);  % HETESSEL
kB_1    =   (fc.^2)      .* kB1_c     +...         % canopy only
    2*fc.*fs     .* kB1_m    +...                  % mixed canopy and soil, 
    (fs.^ 2)     .* kB1_s;                         % soil only

z0h              =   z0m ./ exp(kB_1);                      % roughness height for heat (su 2002 eq 8)
z0h(~isnan(z0h)) = max(z0h(~isnan(z0h)),1.0E-10); % has problem at BSN site of Qomolangma
z0h(z0h>hc)      = 0.125*hc(z0h>hc);    % 
z0h              = abs(z0h); % remove complex data
z0m              = abs(z0m); % remove complex data
d0               = abs(d0); % remove complex data
% return;
%%----------------------- Urban zom
z0m(LCT==13)     = 0.5; % MM5 
z0h(LCT==13)     = 0.5; % MM5
d0(LCT==13)      = 0;   % MM5
%%_______________________ snow and ice
z0m(LCT==15)    = 0.0001; % 
% [z0h(LCT==15) ] = yang_kb_1(z0m(LCT==15),10,10,Uref(LCT==15),LST_K(LCT==15),Tref_K(LCT==15),qaref(LCT==15), Pref(LCT==15));  % Yang et al (2002)
d0(LCT==15)     = 0;
z0h(LCT==15)    = z0m(LCT==15);
%_______________________ water surface
Nu            = 1.327e-5 * (101325 ./ Pref) .* ((Tref_K/ 273.15).^1.81);  % Kinematic viscosity of air (Massman 1999b) (10 cm^2/s)
% compute initial neutral scaling coefficients
[z0m1,ustar1] =  cdntc(Uref,10,Tref_K-273.15);
z0m(LCT==0)   =  z0m1((LCT==0)); % mm5 ,z0m=0.0001;z0h=0.0001;
Re            =  ustar1.*z0m1./Nu; % roughness Reynolds number
z0h(LCT==0)   =  z0m1(LCT==0).*exp(-2.67*Re(LCT==0).^0.25+2.57);
d0(LCT==0)    =  0;
%----------------------  MM5 Perst. Wetland, z0m=0.2;z0h=0.2;
z0m(LCT==11)  =  0.2;
z0h(LCT==11)  =  0.2;
%----------------------  MM5, Permanent ice z0m=0.001,z0h=0.001;
z0m(LCT==23)    = 0.001;

%----------------------- MM5 desert z0m=0.1;z0h=0.1; only valid for Sahra desert
z0m(LCT==26)  =  0.1;
z0h(LCT==26)  =  0.1;
d0(LCT==26)   =  0;

% adjust savanan z0h,
% z0h(LCT==9)  =  z0m(LCT==9);
I=d0+z0m-hc>0;
d0(I)=NaN;
z0m(I)=NaN;
z0h(I)=NaN;

return;