function [z0h]=z0mz0h(zh,nu,ustr,tstr)
%    ##################################################################
%   ##################################################################
%   ######                                                      ######
%   ######                    SUBROUTINE z0mz0h                 ######
%   ######                                                      ######
%   ######                     Developed by                     ######
%   ######     River and Environmental Engineering Laboratory   ######
%   ######                University of Tokyo                   ######
%   ######                                                      ######
%   ##################################################################
%   ##################################################################

%#######################################################################
%   PURPOSE:
%   Surface flux parameterization. Output ustr,tstr,ra,Hsfc
%#######################################################################
%   input
%      real zh        ! reference level of air temperature (m)
%      real z0m       ! aerodynamic roughness length
%      real ustr      ! frictional velocity
%      real tstr      ! =-H/(rhoair*cp*ustr)
%      real nu        ! kinematic viscousity 
%   output
%      real z0h	   ! thermal roughness length

      z0h = 70*nu./ustr.*exp(-7.2*sqrt(ustr).*sqrt(sqrt(abs(-tstr))));  %阳参数化法
%       z0h = min(zh/10,max(z0h,1.0E-5)); % 2015, 1.0E-10 is revised to be 1.0E-5
%       figure;plot(nu ./ ustr,'.');
%       figure;plot(exp(-7.2*sqrt(ustr).*sqrt(sqrt(abs(-tstr)))),'.');