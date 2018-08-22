function [y] = linspace1(d1, d2, n)
%-----------------------------------------------------------------------
%LINSPACE Linearly spaced vector.
%   LINSPACE(x1, x2) generates a row vector of n linearly
%   equally spaced points between x1 and x2.
%
%   See also LOGSPACE, :.
%
%   Copyright 1984-2001 The MathWorks, Inc. 
%   $Revision: 5.11 $  $Date: 2001/04/15 12:02:30 $
%-----------------------------------------------------------------------
[s,t]= size(d1);
d2d1 = ones(s,n-1,'single');
ds   = ones(s,n-1,'single');
d11  = ones(s,n-1,'single');
if length(d1)==1
   y = [d1+(d2-d1)*(0:n-2)/(n-1) d2];
else
%-----------------------方法1    
%     tp(1:n-1)=1; 
%     y = [d1*tp+(d2-d1)*(0:n-2)/(n-1) d2];  
%-----------------------方法2 
    d2md1= d2-d1;
    d2d1= repmat(d2md1,[1,n-1]);% 转成n-1
    n_1 = (0:n-2)/(n-1);
    ds  = repmat(n_1,[s,1]);
    d11 = repmat(d1, [1,n-1]);    % 转成n-1
    y   = d11+d2d1.*ds;
    y   = [y d2];
end