% BOOTSTRAP_MVAR_SPECTRA      Wrapper to bootstrap MVAR_SPECTRAL
% 
%     [out,mvar_boot] = bootstrap_mvar_spectral(v,n,pmin,pmax,nfft,fs,nboot,mvar_point,varargin);
%
%     INPUTS
%     v          - Txn array of signals
%     pmin       - minimum AR order to estimate
%     pmax       - maximum AR order to estimate
%     nfft       - length of FFT
%     fs         - sampling frequency
%     nboot      - # of boostraps
%     mvar_point - point estimate of MVAR statistics
%
%     OPTIONAL
%     varargin   - cell array that gets passed directly to MVAR_SPECTRAL
%
%     OUTPUTS
%     out        - struct array of outputs, depends on analyses requested
%     mvar_boot  - struct array of boostrapped MVAR statistics
%
%     EXAMPLES
%     fgranger_inference_demo

%     $ Copyright (C) 2014 Brian Lau http://www.subcortex.net/ $
%     Released under the BSD license. The license and most recent version
%     of the code can be found on GitHub:
%     https://github.com/brian-lau/Granger
%
%     REVISION HISTORY:
%     brian 11.14.07 written 

function [out,mvar_boot] = bootstrap_mvar_spectral(v,pmin,pmax,nfft,fs,nboot,mvar_point,varargin)

[len,n] = size(v);

for i = 1:nboot
   vboot = v;
   for j = 1:n
      vboot(:,j) = phase_shuffle(v(:,j));
   end
   
   [w,A_hat,Sigma_eps_hat,SBC,FPE,th] = arfit(vboot,pmin,pmax);
   [siglev,res] = arres([0 0]',A_hat,v);
   ind = find(SBC==min(SBC));
   p = pmin:pmax;
   p = p(ind);
   A3D = reshape(A_hat,n,n,p);

   mvar_boot(i) = mvar_spectral(A3D,Sigma_eps_hat,nfft,fs,varargin{:});
end

if sum(strcmp(mvar_boot(1).analyses,'granger'))
   temp = cat(2,mvar_boot.G_yz);
   out.pval_yz = sum(temp>repmat(mvar_point.G_yz,1,nboot),2)/nboot;
   out.pval_yz = out.pval_yz(1:nfft/2);
   temp = cat(2,mvar_boot.G_zy);
   out.pval_zy = sum(temp>repmat(mvar_point.G_zy,1,nboot),2)/nboot;
   out.pval_zy = out.pval_zy(1:nfft/2);
else
   out = [];
end