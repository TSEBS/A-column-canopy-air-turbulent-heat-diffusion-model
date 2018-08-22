%--------------------------------------------------------------------------
%
%                        Canopy_Flow_Model_New_LAI.m
%
%  [RMRS Massman Version 3 = Variable Foliage Distribution -- March 2015]   
%                      
%--------------------------------------------------------------------------
% INPUT  ------- z0h [z0/h]
% z0 = ground surface roughness length and h = canopy height
%__________________________________________________________________________
function [doh,z0oh]=Canopy_Flow_Model(LAI,LCT)

z0h = 0.0025;
%__________________________________________________________________________
% Define z/h vector from top to bottom of the canopy: [0 <= z01 <= 1]
%__________________________________________________________________________
delh = 0.0001;
z01  = 0:delh:1;
nz   = length(z01);
%__________________________________________________________________________
% Define logarithmic profile for bottom portion of the canopy wind profile
%__________________________________________________________________________
uz0  = log(z01/z0h);
%__________________________________________________________________________
% Constants for calculating u*/u(h) = usuh 
%__________________________________________________________________________
vkar = 0.4;
c1   = 0.38;
c3   = 15;
c2   = c1+vkar/log(z0h);
%__________________________________________________________________________
% INPUT ------- Foliage Distribution 
%__________________________________________________________________________
% MassfoldisLAI
% NormfoldisLAI
[hazn]=NormfoldisLAI(LCT); % eq 2 in Massman 2018
% TrifoldisLAI
% UnifoldisLAI
% TrifoldisAS
%__________________________________________________________________________
% CALCULATE -------- dohM(NLAI) and z0ohM(NLAI) 
%__________________________________________________________________________
haz       = LAI*hazn/trapz(hazn); % trapz 梯形数值积分, Eq.1 
%__________________________________________________________________________
% INPUT ------- Drag Coefficient Distribution  
%__________________________________________________________________________
% Massdragdis
[zetaz,hacpz,hacpn]=Massdragdis(haz);
ZETA = hacpn(nz);
%__________________________________________________________________________
% RETURN ----- Normalized Wind Speed Profile ucan [u(z)/u(h)] and 
%__________________________________________________________________________
% canopywind = 'exponent';
% canopywind = 'hypersin';
canopywind = 'hypercos';
[ucan,usuh,nzet,nexp]=MassUProfileLAI(hacpn,nz,zetaz,canopywind);
%__________________________________________________________________________
% COMPUTE ---- Composite wind profile [ubc = u(z)/u(h)] and associated  
%              Composite stress profile [usb = u*(z)/u*(h)] and 
%              d/h and z0/h for whole stand 
%__________________________________________________________________________
inflect = 0;
if inflect == 1
   [ubc,usb,I1,doh,z0oh,rough,lroug] = ... 
               calUUS(ucan,uz0,z01,z0h,delh,vkar,usuh,nz,hacpz);
else      
%  cwd = 'exponent';
   cwd = 'hypercos';
   [ubc,usb,I1,doh,z0oh,rough,lroug] = ... 
               calUUSLAII(ucan,uz0,z01,z0h,delh,vkar,usuh,nz,nzet,nexp,cwd);                    
end

return;
%% __________________________________________________________________________
% RETURN ----- Massman Canopy Wind Speed Ratio = 
%                                  u_canopy(z_0h < z_c/h < 1)/U(H/h+1) 
%
% H     = height above the canopy top in meters 
% h     = tree height (at canopy top) in meters
% Int1  = 1-d/h  
% zfire = desired input height of fire in meters (z0h < zfire/h < 1) 
% 
% UCUH = reduction factor = ratio of wind speed at input height      
%__________________________________________________________________________

h      = 10;
Hflame = 2;
H      = 6.096;
% [Iz0(NLAI) UCUH(NLAI) UCUH2(NLAI) UAUH(NLAI) UAUH2(NLAI)] = ... 
%                       calUCUH(I1,z0oh(NLAI),h,Hflame,H,ubc,delh);
[Iz0(NLAI) UCUH(NLAI) UCUH2(NLAI) UAUH(NLAI) UAUH2(NLAI)] = ... 
                      calUCUHz(I1,z0oh(NLAI),rough,h,Hflame,H,ubc,delh);                  
                    
hn      = 10;
Hflamen = 14;
% [UFUH(NLAI) UFUH2(NLAI) ULUH(NLAI) ULUH2(NLAI)] = ... 
%                         calUFUH(I1,z0oh(NLAI),hn,Hflamen,H); 
[UFUH(NLAI) UFUH2(NLAI) ULUH(NLAI) ULUH2(NLAI)] = ... 
                        calUFUHz(I1,z0oh(NLAI),rough,hn,Hflamen,H); 

