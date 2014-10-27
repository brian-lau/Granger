%     $ Copyright (C) 2014 Brian Lau http://www.subcortex.net/ $
%     Released under the BSD license. The license and most recent version
%     of the code can be found on GitHub:
%     https://github.com/brian-lau/Granger
%
function r = build_r(omega,beta)

p = length(beta);
ind = [1:p];

if max(size(omega)) == 1
   if (omega == 0) || (omega == pi)
      R = [cos(ind.*omega)];
   else
      R = [cos(ind.*omega) ; sin(ind.*omega)];
   end

   r = R*beta;
else
   for i = 1:length(omega)
      r{i} = build_r(omega(i),beta);
   end
   return;
end
