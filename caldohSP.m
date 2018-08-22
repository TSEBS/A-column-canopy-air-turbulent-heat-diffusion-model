%--------------------------------------------------------------------------
function [ Int1,doh,z0oh,rough,lroug ] = caldohSP( delh,vkar,z01,sstress,usuh )
%--------------------------------------------------------------------------
% Use Shaw+Pereira (1982) to calculate the displacement height = d/h = doh 
%     and the surface roughness length = z0oh = z_0/h  
% To determine d/h it is necessary to integrate the momentum flux profile
%     sstress over the canopy depth  
% d/h and z_0/h are used to compute the wind speed ratios 
%__________________________________________________________________________
% rough = 1.07;
rough = 1.25;   
lroug = log(rough);

denom = trapz(sstress);% eq. 15
numer = trapz(z01.*sstress); % eq. 15

doh   = (numer/denom)*(1-sstress(1)); % eq.15
Int1  = 1-doh;
z0oh  = rough*(1-doh)*exp(-vkar/usuh); % eq.16

end