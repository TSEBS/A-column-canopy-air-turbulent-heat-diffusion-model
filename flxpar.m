function [c_u,c_pt]=flxpar(zm,zh,z0m,z0h,wspd,ptsfc,pt1)
% zm=zm2;zh=zh2;ptsfc=pt1;pt1=pt2;
%     PURPOSE:  Compute c_u, c_pt at land surface.
%     The quantity c_u, c_pt are used to obtain surface fluxes for both 
%     the unstable and stable cases.
[lmo]      =  MOlength(zm,zh,z0m,z0h,wspd,ptsfc,pt1);
%以上变量在fortran中统统为六个
[c_u,c_pt] =  CuCpt(lmo,z0m,z0h,zm,zh);
  