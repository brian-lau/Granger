dt = 0.001;
fs = 1/dt;
t = [0:dt:2]';

% noise = randn(length(t),1);
% [B,A] = butter(4,0.25);
% x = filter(B,A,noise);
% v = [x , shiftd(x,0,1,1)];

f1 = 14;
f2 = 75;
x = 4*sin(2*pi*f1*t) + 1*sin(2*pi*f2*t)  + 0.0*randn(length(t),1);
v = [x , shiftd(x,0,1,1)];

nfft = 2^nextpow2(length(x)); % length of fft
l = 1;      % size of first set
m = 1;      % size of second set
n = l+m;    % total size of system
pmin = 1;
pmax = 20;

[w,A_hat,Sigma_eps_hat,SBC,FPE,th] = arfit(v,pmin,pmax);
A3D = reshape(A_hat,n,n,find(SBC==min(SBC)));
mvar_point = mvar_spectral(A3D,Sigma_eps_hat,nfft,fs,'coherence',[],'granger',{l m});

%plot(mvar_point.f,abs(squeeze(mvar_point.S(1,1,:))))

plot(fftshift(mvar_point.f),abs(fft(x,nfft)));
