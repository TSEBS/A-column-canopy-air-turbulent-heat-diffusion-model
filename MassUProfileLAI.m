function [ucan,usuh,nzet,nexp]=MassUProfileLAI(hacpn,nz,zetaz,canopywind)
%__________________________________________________________________________
% Calculate u*/u(h) = usuh 
% Calculate the profile extinction coefficient = nexp
%__________________________________________________________________________
z0h = 0.0025;
vkar = 0.4;
c1   = 0.38;
c3   = 15;
c2   = c1+vkar/log(z0h);
usuh  = c1-c2*exp(-c3*hacpn(nz)); % eq. 10 u*/u(h)
Csurf = 2*usuh*usuh; % eq. 9
nexp  = hacpn(nz)/Csurf; % eq 8
%__________________________________________________________________________
%
% Calculate vertical profiles of wind speed = ucan 
% Choose cosh, sinh, or exponential part of total 
%__________________________________________________________________________
nzet  = nexp*zetaz;
switch canopywind
  case 'hypercos'
    ucan = cosh(nzet)/cosh(nexp); % eq. 12
  case 'hypersin'
    ucan = sinh(nzet)/sinh(nexp);
  case 'exponent'
    ucan = exp(-(nexp-nzet));
end     