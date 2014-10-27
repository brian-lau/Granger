% FGRANGER                   Frequency-domain Granger causality
% 
%     [G_yz,G_zy] = fgranger(H,S,Sigma_eps,l,m);
%
%     INPUTS
%     H         - Transfer matrix, nxnxf
%     S         - Spectral matrix, nxnxf
%     Sigma_eps - Noise covariance, nxn
%     l         - Size of first subsystem (Y)
%     m         - Size of second subsystem (Z), L & M can be excluded if both=1
%
%     OUTPUTS
%     G_yz      - Linear causal influence from Y<-Z
%     G_zy      - Linear causal influence from Z<-Y
%
%     EXAMPLES
%     fgranger_demo

%     $ Copyright (C) 2014 Brian Lau http://www.subcortex.net/ $
%     Released under the BSD license. The license and most recent version
%     of the code can be found on GitHub:
%     https://github.com/brian-lau/Granger
%
%     REVISION HISTORY:
%     brian 08.20.07 written 

function [G_yz,G_zy] = fgranger(H,S,Sigma_eps,l,m)

[mH,nH,nfft] = size(H);

if nargin < 4
   if (mH~=1) && (nH~=1)
      error('H is not 1x1, L and M are required inputs to FGRANGER!');
   else   
      l = 1;
      m = 1;
   end
end

n = l + m;
if (mH~=n) && (nH~=n)
   error('H must be square along the first 2 dimensions in FGRANGER!');
end

lind = 1:l;
mind = (l+1):(l+m);

% Partition the noise covariance
Sigma_eta = Sigma_eps(lind,lind);
Sigma_nu = Sigma_eps(mind,mind);
C = Sigma_eps(lind,mind);

% Linear influence from Z to Y
H_tilde = H;
T = C'*inv(Sigma_eta);
for f = 1:nfft
   H_tilde(lind,lind,f) = H(lind,lind,f) + H(lind,mind,f)*T;
   H_tilde(mind,lind,f) = H(mind,lind,f) + H(mind,mind,f)*T;
end

G_yz = zeros(nfft,1);
for f = 1:nfft
   G_yz(f) = log( det(S(lind,lind,f)) / det(H_tilde(lind,lind,f)*Sigma_eta*H_tilde(lind,lind,f)') );
end

% Linear influence from Y to Z
H_hat = H;
U = C*inv(Sigma_nu);
for f = 1:nfft
   H_hat(lind,mind,f) = H(lind,mind,f) + H(lind,lind,f)*U;
   H_hat(mind,mind,f) = H(mind,mind,f) + H(mind,lind,f)*U;
end

G_zy = zeros(nfft,1);
for f = 1:nfft
   G_zy(f) = log( det(S(mind,mind,f)) / det(H_hat(mind,mind,f)*Sigma_nu*H_hat(mind,mind,f)') );
end

G_yz = real(G_yz);
G_zy = real(G_zy);
