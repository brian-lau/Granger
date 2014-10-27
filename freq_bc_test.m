% FREQ_BC_TEST               Breitung-Candelon frequency-domain test for Granger causality
% 
%     out = freq_bc_test(omega,r,A,v,res,l,m,p,alpha,fs);
%
%     INPUTS
%     omega - Frequencies at which to conduct BC-test. If FS is not provided, this must be specified
%             in *normalized angular frequency*, otherwise it must be specified in angular frequency
%     r     - Null hypothesis values corresponding to tests of the COS and SIN components at
%             frequency OMEGA. Typically testing for 0 amplitude, so r = [0 0]'
%             If r is a 2x2xp array of coefficients, the test will be against the null hypothesis
%             that the amplitude at frequency omega of the model in A is the same as the model in r.
%     A     - nxnxp array of parameters
%     v     - times series, formatted Txn
%     res   - residuals corresponding to AR fit for V
%     l     - dimension of subsystem Y, currently must be 1
%     m     - dimension of subsystem m, currently must be 1
%     p     - AR order
%
%     OPTIONAL
%     alpha - critical level, defaults to 0.05
%     fs    - sampling frequency
%
%     OUTPUTS
%     out   - structure with the fields:
%             .omega_test - frequencies tested (in original units)
%             .pval_yz    - p-values for Y<-Z
%             .F_yz       - F-statistics for Y<-Z
%             .pval_zy    - p-values for Z<-Y
%             .F_zy       - F-statistics for Z<-Y
%
%     EXAMPLES
%     fgranger_demo
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
%     brian 11.21.07 fixed error in passing incorrect SSE for 2<-1
%     brian 11.27.07 added FS argument

function out = freq_bc_test(omega,r,A,v,res,l,m,p,alpha,fs)

if nargin < 10
   fs = 1;
end

if nargin < 9
   alpha = 0.05;
end

if ndims(r) == 3
   A2 = r;
   r = [];
end

% Test for causality at frequency OMEGA_TEST from 1<-2
beta = squeeze(A(1,2,:));
if isempty(r)
   r = build_r(omega,squeeze(A2(1,2,:)));
end
X = format_VAR_design(v(:,2),p);
[pval_yz,F_yz] = bc_test(r,omega,X,beta,sum(res(:,1).^2),alpha,fs);

% Test for causality at frequency OMEGA_TEST from 2<-1
beta = squeeze(A(2,1,:));
if isempty(r)
   r = build_r(omega,squeeze(A2(2,1,:)));
end
X = format_VAR_design(v(:,1),p);
[pval_zy,F_zy] = bc_test(r,omega,X,beta,sum(res(:,2).^2),alpha,fs);

out.omega_test = omega(:);
out.pval_yz = pval_yz(:);
out.F_yz = F_yz(:);
out.pval_zy = pval_zy(:);
out.F_zy = F_zy(:);
