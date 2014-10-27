% BC_TEST                    Breitung-Candelon frequency-domain test for Granger causality
% 
%     [pval,F,Fcrit] = bc_test(r,omega,X,beta,SSE,alpha,fs);
%
%     INPUTS
%     r     - Null hypothesis values corresponding to tests of the COS and SIN components at
%             frequency OMEGA. Typically testing for 0 amplitude, so r = [0 0]'
%     omega - Frequencies at which to conduct BC-test. If FS is not provided, this must be specified
%             in *normalized angular frequency*, otherwise it must be specified in angular frequency
%     X     - Portion of the AR design matrix corresponding to the signal one is testing the 
%             causality from
%     beta  - Portion of the coefficient vector corresponding to the signal one is testing the 
%             causality from
%     SSE   - sum of squared errors of full model
%
%     OPTIONAL
%     alpha - critical level, defaults to 0.05
%     fs    - sampling frequency
%
%     OUTPUTS
%     pval  - corresponding to tested frequencies
%     F     - corresponding to tested frequencies
%     Fcrit - corresponding to ALPHA
%
%     SEE ALSO
%     freq_bc_test, which is just a handy wrapper for calling bc_test
%
%     REFERENCE
%     Breitung & Candelon (2006). Testing for short- and long-run causality: A frequency-domain 
%     approach. Journal of Econometrics, 132(2), 363-378. 
%     http://dx.doi.org/10.1016/j.jeconom.2005.02.004

%     $ Copyright (C) 2014 Brian Lau http://www.subcortex.net/ $
%     Released under the BSD license. The license and most recent version
%     of the code can be found on GitHub:
%     https://github.com/brian-lau/Granger
%
%     REVISION HISTORY:
%     brian 11.14.07 written 
%     brian 11.21.07 cleaned up handling on 0 and pi
%     brian 11.27.07 added FS argument
%     brian 11.01.11 avoid inverting X'X repeatedly 

function [pval,F,Fcrit] = bc_test(r,omega,X,beta,SSE,alpha,fs,invFlag);

if nargin >= 7
   if ~isempty(fs)
      % Convert standard frequency (Hz) to normalized angular frequency
      omega = (2*pi)*(omega./fs);
   end
end

if nargin < 6
   alpha = 0.05;
end
if nargin < 8
   invFlag = false;
end

p = length(beta);
T = size(X,1);

if max(size(omega)) == 1
   ind = [1:p];
   
   if (omega == 0) | (omega == pi)
      R = [cos(ind.*omega)];
      r = r(1);
   else
      R = [cos(ind.*omega) ; sin(ind.*omega)];
   end
   if invFlag
      % inv(X'*X) passed in as X
      F = (r - R*beta)'*inv(R*X*R')*(r - R*beta);
      T = invFlag;
   else
      F = (r - R*beta)'*inv(R*inv(X'*X)*R')*(r - R*beta);
   end
   q = rank(R);
   F = (F/q) / (SSE/(T-q*p));
   pval = 1 - fcdf(F,q,T-q*p);
   
   if nargout == 3
      Fcrit = finv(1-alpha/2,q,T-q*p);
   end
else
   invXX = inv(X'*X);
   for i = 1:length(omega)
      if ~iscell(r)
         [pval(i),F(i)] = bc_test(r,omega(i),invXX,beta,SSE,alpha,[],size(X,1));
      else
         [pval(i),F(i)] = bc_test(r{i},omega(i),invXX,beta,SSE,alpha,[],size(X,1));
      end
   end
   return
end