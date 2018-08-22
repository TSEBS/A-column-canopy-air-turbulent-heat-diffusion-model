function [psim] = PSIm(zeta)
Y                                   =   -zeta;                  % Brutsaert 2008, p50
% Integrated stability function 1a. Stability correction function for momentum, eq.(16)
% Y=-z/L, z is the height, L the Obukhov length, both specified in the calling statement.
% For stable conditions we use the expressions proposed by Beljaars and Holtslag (1991)and evaluated by Van den Hurk and
% Holtslag (1995) can be used.
a_s     =   1.0;                                                % constants, p. 122, van der Hurk & Holtslag (1997)
b_s     =   0.667;                                              % constants, p. 122, van der Hurk & Holtslag (1997)
c_s     =   5.0;                                                % constants, p. 122, van der Hurk & Holtslag (1997)
d_s     =   0.35;                                               % constants, p. 122, van der Hurk & Holtslag (1997)% QUESTION: In page 24, d_s=1 
a_u     =   0.33;                                               % constants, p. 49, Brutsaert(2008)
b_u     =   0.41;                                               % constants, p. 49, Brutsaert(2008)
I_s     =   (Y < 0.0);                                          % STABLE conditions (According to Beljaars & Holtslag, 1991, eq. 13)
I_u     =   (Y >= 0.0);                                         % UNSTABLE conditions(% According to Brutsaert 2008, p50)
I_u1    =   (Y <= b_u^(-3));
I_u2    =   (Y >  b_u^(-3));
y_s     =   -Y;                                                 % due to formulation of Beljaars and Holtslag 1991
y_u1    =   Y;
y_u2    =   b_u^(-3);
x_u1    =   (y_u1/a_u).^(1/3);    
x_u2    =   (y_u2/a_u).^(1/3);        
y_u     =   I_u1.*y_u1 + I_u2.*y_u2;
x_u     =   I_u1.*x_u1 + I_u2.*x_u2;   
PSI0    =   -log(a_u) + sqrt(3)*b_u*(a_u^(1/3))*pi/6; 
psim_s  =   -(a_s * y_s + b_s*(y_s - c_s/d_s).*exp(-d_s*y_s)+b_s*c_s/d_s);
tp1     =   b_u*a_u^(1/3)/2 * log((1+x_u).^2 ./ (1-x_u+x_u.^2));
tp2     =   sqrt(3) * b_u*a_u^(1/3) * atan((2*x_u - 1)/sqrt(3));
psim_u  =   log(a_u+y_u) - 3*b_u*(y_u.^(1/3)) +  tp1 + tp2 + PSI0;
psim    =   I_s .* psim_s + I_u .* psim_u;
return