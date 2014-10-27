function [pval,F] = test_freq_test(omega_test);

omega = pi/2;

alpha = 0.3;
A1 = [ 0.0   alpha;
       0.0   0.0 ];
A2 = [ 0.0  -2*alpha*cos(omega);
       0.0  -0.0 ];
A3 = [ 0.0  alpha;
       0.0  0.0 ];
Sigma_eps = [ 1.00  0.00;
              0.00  1.00 ];

% No causality
% alpha = 0;
% A1 = [ 0.1   alpha;
%        0.0   1 ];
% A2 = [ 0.0  -2*alpha*cos(omega);
%        0.0  -0.1 ];
% A3 = [ 0.0  alpha;
%        0.0  0.0 ];
% Sigma_eps = [ 0.50  0.20;
%               0.20  0.50 ];

% % This model is finicky?
% alpha = 0.3;
% A1 = [ 0.1   alpha;
%        -1.0   0.1 ];
% A2 = [ 0.0  -2*alpha*cos(omega);
%        0.0  -0.2 ];
% A3 = [ 0.0  alpha;
%        0.0  0.3 ];
% Sigma_eps = [ 0.5  0.2;
%               0.2  0.5 ];

A = [ A1 A2 A3 ];

len = 300;
v = arsim([0 0],A,Sigma_eps,len);

l = 1;      % size of first set
m = 1;      % size of second set
n = l+m;    % total size of system
pmin = 3;
pmax = 3;
nfft = 2^nextpow2(len); % length of fft
fs = 2*pi;

[w,A_hat,Sigma_eps_hat,SBC,FPE,th] = arfit(v,pmin,pmax);
%[w,A_hat,Sigma_eps_hat,SBC,AIC] = arfit2(v,pmin,pmax);
[siglev,res] = arres([0 0]',A_hat,v);
ind = find(SBC==min(SBC));
p = pmin:pmax;
p = p(ind);
A3D = reshape(A_hat,n,n,p);
%A3D = reshape(A,n,n,p);

mvar_point = mvar_spectral(A3D,Sigma_eps_hat,nfft,fs,'coherence',[],'granger',{l m});
figure; hold on
plot(fftshift(mvar_point.f),real(mvar_point.G_yz),'r')
plot(fftshift(mvar_point.f),real(mvar_point.G_zy),'k--')

% ind = [1:p];
% R = [cos(ind.*omega) ; sin(ind.*omega)];
%beta = squeeze(A3D(1,2,:));
beta = squeeze(A3D(2,1,:));
r = [0 0]';
%keyboard
X = format_VAR_design(v(:,1),p);
% F = (r - R*beta)'*inv(R*inv(X'*X)*R')*(r-R*beta)
% 
% pval = 1 - fcdf(F,2,size(X,1)-2*p)

[pval,F] = bc_test(r,omega_test,X,beta,sum(res(:,1).^2));
