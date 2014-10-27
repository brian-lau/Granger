% MINOR                      Matrix minor
% 
%     m = minor(M,i,j);
%
%     INPUTS
%     M - matrix
%     i - row index
%     j - column index
%
%     OUTPUTS
%     m - determinant
%

%     $ Copyright (C) 2014 Brian Lau http://www.subcortex.net/ $
%     Released under the BSD license. The license and most recent version
%     of the code can be found on GitHub:
%     https://github.com/brian-lau/Granger
%
%     REVISION HISTORY:
%     brian 08.20.07 written 

function m = minor(M,i,j);

[m,n] = size(M);

if m ~= n
   error('Input must be square for MINOR');
end

M(i,:) = [];
M(:,j) = [];
m = det(M);
