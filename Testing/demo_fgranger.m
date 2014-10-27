% Bivariate system, with causality in only one direction
%%%%%%%%%%%%%

l = 1;      % size of first set
m = 1;      % size of second set
n = l+m;    % total size of system
p = 3;      % order of MVAR
fs = 1000;  % sampling frequency
nfft = fs;  % length of fft
len = 1000; % number of samples

omega = pi/2;
beta = 0.3;
% Coefficients for tau=1
A1 = [ 0.1   beta;
       0   0.1 ];
% Coefficients for tau=2
A2 = [ 0.0  -2*beta*cos(omega);
       0.0  -0.2 ];
% Coefficients for tau=3
A3 = [ 0.0  beta;
       0.0  0.3 ];
% Noise covariance
Sigma_eps = [ 0.50  0.20;
              0.20  0.50 ];
% Concatonate parameters
A = [ A1 A2 A3 ];

% Simulate using function from ARFIT package
v = arsim([0 0],A,Sigma_eps,len);

% Fit AR model and reshape coefficient matrix, MVAR_SPECTRAL expects the coefficient matrix to be 3-D, nxnxp
[w,A_hat,Sigma_eps_hat,SBC,FPE,th] = arfit(v,p,p);
[siglev,res] = arres([0 0]',A_hat,v);
A3D_hat = reshape(A_hat,n,n,p);
%A3D = reshape(A,n,n,p);

mvar = mvar_spectral(A3D_hat,Sigma_eps_hat,nfft,fs,'coherence',[],'granger',{l m});

figure;
subplot(2,2,1);
plot(mvar.f(1:fs/2),abs(squeeze(mvar.S(1,1,1:fs/2))),'r');
title('Spectral power of signal 1')
subplot(2,2,2); hold on
plot(mvar.f(1:fs/2),abs(squeeze(mvar.R(1,2,1:fs/2))),'r');
plot(mvar.f(1:fs/2),mvar.G_yz(1:fs/2),'b');
title({'Coherence of signal 1, red' 'Granger causality (1<-2), blue'})
subplot(2,2,3); hold on
plot(mvar.f(1:fs/2),abs(squeeze(mvar.R(1,2,1:fs/2))),'r');
plot(mvar.f(1:fs/2),mvar.G_zy(1:fs/2),'b');
title({'Coherence of signal 2, red' 'Granger causality (2<-1), blue'})
subplot(2,2,4);
plot(mvar.f(1:fs/2),abs(squeeze(mvar.S(2,2,1:fs/2))),'r');
title('Spectral power of signal 2')

%----- Hypothesis testing
%--- Time-domain
% Basic test against null that all parameters are zero, in this case
% it is easy to set the contrast up directly
H = eye(n*n*p);
h = zeros(n*n*p,1);
pval = mvar_tgranger_test(v,A3D_hat,H,h)

% Test against null that there is no causality from 2<-1
clear R r;
R{1} = zeros(n,n,p);
R{1}(2,1,1) = 1;
r{1} = 0;
R{2} = zeros(n,n,p);
R{2}(2,1,2) = 1;
r{2} = 0;
pval = mvar_tgranger_test(v,A3D_hat,R,r)

% Test against null that there is no causality from 1<-2
clear R r;
R{1} = zeros(n,n,p);
R{1}(1,2,1) = 1;
r{1} = 0;
R{2} = zeros(n,n,p);
R{2}(1,2,2) = 1;
r{2} = 0;
pval = mvar_tgranger_test(v,A3D_hat,R,r)

%--- Frequency-domain
% Using a test developed by Breitung & Candelon (2006)
alpha = 0.05;
bc = freq_bc_test(mvar.f(1:fs/2),[0 0]',A3D_hat,v,res,l,m,p,alpha,fs);

figure;
subplot(2,2,2); hold on
title({'Granger causality (1<-2), blue' 'F-statistic (1<-2), green' 'Significant causality marked with x''s'})
plotyy(mvar.f(1:fs/2),mvar.G_yz(1:fs/2),mvar.f(1:fs/2),bc.F_yz);
temp = mvar.G_yz(1:fs/2);
temp(bc.pval_yz>=alpha) = NaN;
plot(mvar.f(1:fs/2),temp,'rx')
subplot(2,2,3); hold on
title({'Granger causality (2<-1), blue' 'F-statistic (2<-1), green' 'Significant causality marked with x''s'})
plotyy(mvar.f(1:fs/2),mvar.G_zy(1:fs/2),mvar.f(1:fs/2),bc.F_zy)
temp = mvar.G_zy(1:fs/2);
temp(bc.pval_zy>=alpha) = NaN;
plot(mvar.f(1:fs/2),temp,'rx')

% Test using a brute-force bootstrap (fourier-shuffled)
nboot = 100;
pmin = 1;
pmax = 10;
tic;
boot = bootstrap_mvar_spectral(v,pmin,pmax,nfft,fs,nboot,mvar,'granger',{l m});
toc

figure;
subplot(2,2,2); hold on
title({'Granger causality (1<-2), blue' 'Significant causality (bootstrap) marked with +''s'})
plot(mvar.f(1:fs/2),mvar.G_yz(1:fs/2));
temp = mvar.G_yz(1:fs/2);
temp(boot.pval_yz>=alpha) = NaN;
plot(mvar.f(1:fs/2),temp,'r+')
subplot(2,2,3); hold on
title({'Granger causality (2<-1), blue' 'Significant causality (bootstrap) marked with +''s'})
plot(mvar.f(1:fs/2),mvar.G_zy(1:fs/2))
temp = mvar.G_zy(1:fs/2);
temp(boot.pval_zy>=alpha) = NaN;
plot(mvar.f(1:fs/2),temp,'r+')

% % Similar example for multivariate system
% %%%%%%%%%%%%%
% A1 = [ 0.4 0.31 0.3 0.3;
%      -.2 0.4 0.0 0.0;
%      0.0 0.0 0.4 0.0;
%      0.0 0.0 0.0 0.4];
% A2 = [ 0.3 0.5 0.1 0.6;
%      -.1 0.3 0.0 0.0;
%      0.0 0.0 0.3 0.0;
%      0.0 0.0 0.0 0.3];
% Sigma_eps = [ 1.0 0.25 0.0 0.0;
%               0.25 1.0 0.0 0.0;
%               0.0 0.0 1.0 0.25;
%               0.0 0.0 0.25 1.0];
% A = [ A1 A2 ];
% 
% v = arsim([0 0 0 0],A,Sigma_eps,200);
% 
% nfft = 256;
% l = 2;
% m = 2;
% n = l+m;
% p = 2;
% 
% fs = 2*pi;
% 
% A3D = reshape(A,n,n,p);
% mvar = mvar_spectral(A3D,Sigma_eps,nfft,fs,'granger',{l m});
% 
% figure;
% subplot(2,2,2); hold on
% plot(mvar.f,mvar.G_yz,'r');
% subplot(2,2,3); hold on
% plot(mvar.f,mvar.G_zy,'r');
