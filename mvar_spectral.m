% MVAR_SPECTRAL              Spectral analyses for multivariate autoregressive models
% 
%     mvar = mvar_spectral(A,Sigma_eps,p,nfft,fs,varargin);
%
%     INPUTS
%     A         - MVAR coefficient matrix, nxnxp
%     Sigma_eps - Noise covariance, nxn
%     nfft      - length of FFT
%     fs        - sampling frequency
%
%     OPTIONAL
%     varargin  - cell array where requested analysis is followed by parameters (or
%                 empty array if none)
%                 available analyses:
%                   'coherence', no options
%                   'partialcoherence', no options
%                   'partialdirectedcoherence', no options
%                   'directedtransferfunction', no options
%                   'fgranger', {l m} indicate the size of the 1st and 2nd subsystems
%
%     OUTPUTS
%     mvar      - struct 
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

function mvar = mvar_spectral(A,Sigma_eps,nfft,fs,varargin);

if rem(length(varargin),2) ~= 0
   error('Bad VARARGIN to MVAR_SPECTRAL!');
end

count = 1;
for i = 1:2:length(varargin)
   analyses{count} = varargin{i};
   params{count} = varargin{i+1};
   count = count + 1;
end

[m,n] = size(Sigma_eps);
if m~=n 
   error('SIGMA_EPS must be square for MVAR_SPECTRAL!')
end

B = cat(3,eye(n),-A);
fB = fft(B,nfft,3);

H = zeros(size(fB));
for f = 1:nfft
   H(:,:,f) = inv(fB(:,:,f));
end

S = zeros(size(H));
for f = 1:nfft
   S(:,:,f) = H(:,:,f)*Sigma_eps*H(:,:,f)';
end

mvar.A = A;
mvar.B = B;
mvar.Sigma_eps = Sigma_eps;
mvar.nfft = nfft;
mvar.f = freq(fs,nfft);
mvar.H = H;
mvar.S = S;
mvar.analyses = analyses;
mvar.params = params;

for analysis = 1:length(analyses)
   if strcmp(analyses{analysis},'coherence')
      R = repmat(eye(n),[1 1 nfft]);
      for i = 1:n
         for j = (i+1):n
            R(i,j,:) = squeeze(S(i,j,:)) ./ sqrt(squeeze(S(i,i,:)).*squeeze(S(j,j,:)));
         end
      end
      mvar.R = R;
   end

   if strcmp(analyses{analysis},'partialcoherence')
      D = repmat(eye(n),[1 1 nfft]);
      for i = 1:n
         for j = (i+1):n
            for f = 1:nfft
               D(i,j,f) = minor(S(:,:,f),i,j) ./ sqrt(minor(S(:,:,f),i,i).*minor(S(:,:,f),j,j));
            end
         end
      end
      mvar.D = D;
      
      % Baccala & Sameshima (2001, eq 14)
      % D = repmat(eye(n),[1 1 nfft]);
      % temp = inv(Sigma_eps);
      % for i = 1:n
      %    for j = (i+1):n
      %       for f = 1:nfft
      %          D(i,j,f) = fB(:,i,f)'*temp*fB(:,j,f) ./ sqrt((fB(:,i,f)'*temp*fB(:,i,f)) * (fB(:,j,f)'*temp*fB(:,j,f)));
      %       end
      %    end
      % end
   end
   
   if strcmp(analyses{analysis},'directedtransferfunction')
      DTF = repmat(eye(n),[1 1 nfft]);
      for i = 1:n
         for j = i:n
            for f = 1:nfft
               DTF(i,j,f) = H(i,j,f) ./ sqrt(H(i,:,f)*H(i,:,f)');
               %DTF(i,j,f) = H(i,j,f) ./ sqrt(sum(abs(H(i,:,f)).^2));
            end
         end
      end
      mvar.DTF = DTF;
   end
   
   if strcmp(analyses{analysis},'partialdirectedcoherence')
      PDC = repmat(eye(n),[1 1 nfft]);
      for i = 1:n
         for j = i:n
            for f = 1:nfft
               PDC(i,j,f) = fB(i,j,f) ./ sqrt(fB(i,:,f)*fB(i,:,f)');
            end
         end
      end
      mvar.PDC = PDC;
   end

   if strcmp(analyses{analysis},'fgranger') || strcmp(analyses{analysis},'granger')
      [G_yz,G_zy] = fgranger(H,S,Sigma_eps,params{analysis}{1},params{analysis}{2});
      mvar.G_yz = G_yz;
      mvar.G_zy = G_zy;
   end
end

return

