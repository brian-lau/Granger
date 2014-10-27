% Bivariate system, with causality in only one direction
len = 1000;
nboot = 10;
%%%%%%%%%%%%%
% Coefficients for tau=1
A1 = [ 0.4   0.7;
       0.0   0.7 ];
% Coefficients for tau=2
A2 = [ 0.35  -0.7;
       0.0  -0.5 ];
% Noise covariance
Sigma_eps_true = [ 1.00  0.10;
              0.10  1.50 ];
A_true = [ A1 A2 ];
% Simulate using function from ARFIT package
v = arsim([0 0],A_true,Sigma_eps,len);


v = randn(len,1);
v = [v , shiftd(v,0,1,1) + 0.0*randn(len,1)];

nfft = len; % length of fft
l = 1;      % size of first set
m = 1;      % size of second set
n = l+m;    % total size of system
fs = 2*pi;  % sampling frequency, chosen so frequencies go -pi:pi
pmin = 1;
pmax = 5;

[w,A_hat,Sigma_eps_hat,SBC,FPE,th] = arfit(v,pmin,pmax);
A3D = reshape(A_hat,n,n,find(SBC==min(SBC)));
mvar_point = mvar_spectral(A3D,Sigma_eps_hat,nfft,fs,'coherence',[],'granger',{l m});
for i = 1:nboot
   if 1
      ind = unidrnd(len,len,n);
   else
      ind = unidrnd(len,len,1);
      ind = [ind , ind];
   end
   [w,A_hat,Sigma_eps_hat,SBC,FPE,th] = arfit(v(ind),pmin,pmax);
   p(i) = find(SBC==min(SBC));      % order of MVAR
   % % MVAR_SPECTRAL expects the coefficient matrix to be 3-D, nxnxp
   A3D = reshape(A_hat,n,n,p(i));
   mvar(i) = mvar_spectral(A3D,Sigma_eps_hat,nfft,fs,'coherence',[],'granger',{l m});
end