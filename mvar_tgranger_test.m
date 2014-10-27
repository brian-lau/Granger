% MVAR_TGRANGER_TEST         Wald test for general linear hypotheses on AR parameters
% 
%     [pval,wald] = mvar_tgranger_test(X,A,R,r,convmtxform);
%
%     INPUTS
%     X           - Txn array of signals
%     A           - nxnxp array of parameters
%     R           - Contrast matrix. Since the parameter matrix A is 3-dimensional, R should
%                   be a cell array where each element is a nxnxp matrix specifying a single
%                   contrast. It is possible to also pass in a matrix, but you must understand
%                   how the elements must be arranged to do this (see code).
%     r           - Null hypothesis values corresponding to contrasts in R. If R is a cell array,
%                   r must also be a cell array, if R is a matrix, r must be a vector
%
%     OPTIONAL
%     convmtxform - Binary flag indicating form of XDOT matrix (mainly for testing),
%                   defaults to 0 (Hamilton form)
%
%     OUTPUTS
%     pval        - P-value
%     wald        - Wald statistic
%
%     EXAMPLES
%     fgranger_demo

%     $ Copyright (C) 2014 Brian Lau http://www.subcortex.net/ $
%     Released under the BSD license. The license and most recent version
%     of the code can be found on GitHub:
%     https://github.com/brian-lau/Granger
%
%     REVISION HISTORY:
%     brian 11.14.07 written 

function [pval,wald] = mvar_tgranger_test(X,A,R,r,convmtxform)

if nargin < 5
   convmtxform = false;
end

p = size(A,3); % VAR order
n = size(A,1); % # of signals
T = length(X); % # of samples/signal

if ~convmtxform
   % Formatting is per Hamilton (1994), page 298
   % Xdot(t,:) = [x1(t-1) x2(t-1) ... xn(t-1) , x1(t-2) x2(t-2) ... xn(t-2) , ... , x1(t-p) x2(t-p) ... xn(t-p)]
   % Lambda(:,i) = [A(i,1,1) A(i,2,2) ... A(i,n,p)]'
   Xdot = format_VAR_design(X,p,convmtxform);
   
   Lambda = zeros(n*p,n);
   for i = 1:n
      temp = A(i,:,:);
      Lambda(:,i) = temp(:);
   end

   if iscell(R) && iscell(r)
      H = zeros(length(r),n*n*p);
      h = zeros(length(r),1);
      temp2 = zeros(n*p,n);
      for i = 1:length(r)
         for j = 1:n
            temp = R{i}(j,:,:);
            temp2(:,j) = temp(:);
         end
         H(i,:) = temp2(:)';
         h(i) = r{i};
      end
   else
      H = R;
      h = r;
   end
else
   % My version, with convolution matrices for each variable stacked
   % Xdot(t,:) = [x1(t-1) x1(t-2) ... x1(t-p) , x2(t-1) x2(t-2) ... x2(t-p) , ... , xn(t-1) xn(t-2) ... xn(t-p)]
   % Lambda(:,i) = [A(i,1,:) A(i,2,:) ... A(i,n,:)]'
   Xdot = format_VAR_design(X,p,convmtxform);

   Lambda = zeros(n*p,n);
   for i = 1:n
      temp = A(i,:,:);
      Lambda(:,i) = reshape(squeeze(temp)',n*p,1);
   end

   if iscell(R) && iscell(r)
      H = zeros(length(r),n*n*p);
      h = zeros(length(r),1);
      temp2 = zeros(n*p,n);
      for i = 1:length(r)
         for j = 1:n
            temp = R{i}(j,:,:);
            temp2(:,j) = reshape(squeeze(temp)',n*p,1);
         end
         H(i,:) = temp2(:)';
         h(i) = r{i};
      end
   else
      H = R;
      h = r;
   end
end

% Residuals
res = X - (Lambda'*Xdot')';
% Asymptotic parameter variances
Sigma_hat = (res'*res)/(T-n*p);
var_vecLambda = kron(Sigma_hat,inv(Xdot'*Xdot));

wald = (H*Lambda(:) - h)'*inv(H*var_vecLambda*H')*(H*Lambda(:) - h);
pval = 1 - chi2cdf(wald,rank(H));

