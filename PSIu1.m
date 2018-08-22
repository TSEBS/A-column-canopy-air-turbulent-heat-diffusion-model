function [psiu] = PSIu1(hc,LAI,Zref,zstar,d0,L,Id)
% Integrated stability function
% 1a. Stability correction function for momentum, eq.(16)
% Y=-z/L, z is the height, L the Obukhov length, both specified in the calling statement.
% For stable conditions
global Cd;global phiuc;
stp      =  200;     %
tp(1:stp)=  1;       %
k        =  0.4;
% xi      = (Zref-d0)./(zstar-d0);
beta0    = 0.3;  % beta0=ustar/uh; uh wind speed at the top of the canopy, beta0=0.3, a typical value under neutral conditions.
% Physick and Garatt (1995),phym =  0.5*exp(0.7*xi);
%    z       =  linspace1(Zref,zstar,stp);  % zstar>Zref.
%    xi      =  (z-d0*tp)./(zstar*tp-d0*tp);
z       =  linspace1(Zref-d0,zstar-d0,stp);  % zstar>Zref.integration goes from a certain height in the roughness sublayer up to the top of the roughness sublayer.¼û Molder 1999, AFM,
xi      =  z./(zstar*tp-d0*tp);% Kitsri thesis
L       =  L*tp;    % transfer from m*1 to m*n matrix
phym    =  0.5*exp(0.7*xi);    % eq 18 Physick and Garratt 1995, BLM,
fz      =  (1-phym)./z;        %  see equation 3-4 in Kitsri thesis,
if Id==1 % for momentum
    fz(L<0) =  (1-15*z(L<0)./L(L<0)).^(-0.25).*(1-phym(L<0))./z(L<0); % unstable,Andreas, E.L., Estimating Cn2 Over Snow And Sea Ice From Meteorological Quantities. 1988
    fz(L>0) =  (1+4.7*z(L>0)./L(L>0)).*(1-phym(L>0))./z(L>0); % stable,Andreas, E.L., Estimating Cn2 Over Snow And Sea Ice From Meteorological Quantities. 1988
    %            psiu(i,j) = double(int(fz(i,j),Zref(1,1),zstar)); %equation 3-4 in Kitsri thesis.integration of phym from z to z*,see equation 3-4 in Kitsri
    psiu    =  sum(fz,2).*(zstar-Zref)/stp;% equation 3-4 in Kitsri thesis.integration of phym from z to z*,see equation 3-4 in Kitsri
else % for heat
    Pr      =  0.71;
    fz(L<0) =  Pr*(1-9*z(L<0)./L(L<0)).^(-0.5).*(1-phym(L<0))./z(L<0); % unstable,Andreas, E.L., Estimating Cn2 Over Snow And Sea Ice From Meteorological Quantities. 1988
    fz(L>0) =  (Pr+4.7*z(L>0)./L(L>0)).*(1-phym(L>0))./z(L>0); % stable,Andreas, E.L., Estimating Cn2 Over Snow And Sea Ice From Meteorological Quantities. 1988
    psiu    =  sum(fz,2).*(zstar-Zref)/stp;% equation 3-19 in Kitsri thesis.integration of phym from z to z*,see equation 3-4 in Kitsri
end