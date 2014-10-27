% FORMAT_VAR_DESIGN          Format VAR design matrix
% 
%     Xdot = format_VAR_design(X,p,convmtxform);
%
%     INPUTS
%     X           - Txn array of signals
%     p           - AR order
%
%     OPTIONAL
%     convmtxform - set to true (1) if you want Xdot in convolution matrix form
%                   (see the code, mostly for testing), defaults to 0, which is
%                   the standard (Hamilton, 1994) format for the design matrix
%
%     OUTPUTS
%     Xdot        - design matrix
%
%     EXAMPLES

%     $ Copyright (C) 2014 Brian Lau http://www.subcortex.net/ $
%     Released under the BSD license. The license and most recent version
%     of the code can be found on GitHub:
%     https://github.com/brian-lau/Granger
%
%     REVISION HISTORY:
%     brian 11.21.07 written 
%     brian 11.28.07 fixed an indexing error for p ~= n

function Xdot = format_VAR_design(X,p,convmtxform)

if nargin < 3
   convmtxform = false;
end

n = size(X,2); % # of signals
T = length(X); % # of samples/signal

if ~convmtxform
   % Formatting is per Hamilton (1994), page 298
   % Xdot(t,:) = [x1(t-1) x2(t-1) ... xn(t-1) , x1(t-2) x2(t-2) ... xn(t-2) , ... , x1(t-p) x2(t-p) ... xn(t-p)]
   Xdot_temp = zeros(T,n*p);
   for i = 1:n
      ind = (i-1)*p + (1:p);
      temp = convmtx(X(:,i),p);
      Xdot_temp(:,ind) = [zeros(1,p) ; temp(1:end-p,:)];
   end
   if n == 1
      Xdot = Xdot_temp;
   else
      Xdot = zeros(T,n*p);
      for i = 1:p
         ind1 = (i-1)*n + (1:n);
         ind2 = i:p:(n*p);
         Xdot(:,ind1) = Xdot_temp(:,ind2);
      end
   end
else
   % My version, with convolution matrices for each variable stacked
   % Xdot(t,:) = [x1(t-1) x1(t-2) ... x1(t-p) , x2(t-1) x2(t-2) ... x2(t-p) , ... , xn(t-1) xn(t-2) ... xn(t-p)]
   % Lambda(:,i) = [A(i,1,:) A(i,2,:) ... A(i,n,:)]'
   Xdot = zeros(T,n*p);
   for i = 1:n
      ind = (i-1)*p + (1:p);
      temp = convmtx(X(:,i),p);
      Xdot(:,ind) = [zeros(1,p) ; temp(1:end-p,:)];
   end
end