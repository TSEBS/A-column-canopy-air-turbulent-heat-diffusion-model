function [psih] = PSIh(zeta) 
Y      =   -zeta;                                                %Brutsaert 2008, p50
% Integrated stability function 1b.
% Stability correction function for heat, eq.(17) Y=-z/L, z is the height, L the Obukhov length, both specified in the calling statement.
% For stable conditions we use the expressions proposed by Beljaars and Holtslag (1991)
% and evaluated by Van den Hurk and Holtslag (1995) can be used.
a_s     =   1.0;                                                % constants, p. 122, van der Hurk & Holtslag (1995)
b_s     =   0.667;                                              % constants, p. 122, van der Hurk & Holtslag (1995)
c_s     =   5.0;                                                % constants, p. 122, van der Hurk & Holtslag (1995)
d_s     =   0.35;                                               % constants, p. 122, van der Hurk & Holtslag (1995)% QUESTION: d_s=1 in page 34 of SEBS document of Su
c_u     =   0.33;                                               % constants, p. 443, Brutsaert, 2008
d_u     =   0.057;                                              % constants, p. 443, Brutsaert, 2008
n       =   0.78;                                               % constants, p. 443, Brutsaert, 2008
y_s     =   -Y;                                                 % due to formulation of Beljaars and Holtslag 1991
y_u     =   Y;
I_s     =   (Y < 0);
I_u     =   (Y >= 0);
% STABLE conditions (According to Beljaars & Holtslag, 1991 eq. 13    )    
tp1     =  (1 + 2 * a_s/3 * y_s).^1.5;
tp2     =  b_s * (y_s - c_s /d_s).*exp(-d_s * y_s);
psih_s  =   -(tp1 + tp2 + b_s * c_s /d_s - 1);
% UNSTABLE conditions (According to Brutsaert 2008, p50    
psih_u  =   ((1 - d_u) / n) * log((c_u + (y_u.^ n)) / c_u);     
psih    =   I_s.*psih_s + I_u.*psih_u;
return