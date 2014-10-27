% TODO
%
% do model selection?
% standardize by std?
% x handle different length trials
%   this will require checking and tossing a trial if data doesn't extend
%   throuh the whole window. Output a count, and there should be a
%   threshold number otherwise output an empty struct

% winLen and stepLen are in SAMPLES to avoid rounding errors

%function [mvar,bc,winStartT] = sliding_mvar_spectral(dat,l,m,winLen,stepLen,flag,pmin,pmax,nfft,fs,demean);
function [mvar,bc,winStartT] = sliding_mvar_spectral(dat,l,m,fs,winLen,stepLen,pmin,pmax,overrides);

if exist('overrides','var')
   for i = 1:length(overrides)
      eval([overrides{i} ';']);
   end
end

if ~exist('fitFlag','var')
   fitFlag = 'arfit2';
end
if ~exist('demean','var')
   demean = true;
end
if ~exist('alpha','var')
   alpha = 0.05;  % critical level for statistical tests
end
if ~exist('nfft','var')
   nfft = fs;
end
if ~exist('omega','var')
   f = freq(fs,nfft);
   omega_vec = f(1:fs/2);
else
   f = freq(fs,nfft);
   ind = (f>=omega(1)) & (f<=omega(2));
   omega_vec = f(ind);
end

n = l + m;     % total size of system
nTrials = length(dat);

for i = 1:nTrials
   nSamples(i) = length(dat(i).v);
end
maxSamples = max(nSamples);

% Get window boundaries
winStart = 1:stepLen:maxSamples;
winEnd = winStart + winLen - 1;
ind = winEnd > maxSamples;
winStart(ind) = [];
winEnd(ind) = [];
nWin = length(winStart);

for i = 1:nWin
   % make a function, and deal with different data lengths per trial
   datChunk = format_arfit_input(dat,winStart(i),winEnd(i),fitFlag,demean);
   %keyboard
   %size(datChunk)
   if ~isempty(datChunk) % probably need to check that it is bigger than pmin
      % Fit AR model and reshape coefficient matrix
      switch fitFlag
         case 'arfit'
            [w,A_hat,Sigma_eps_hat,SBC] = arfit(datChunk,pmin,pmax);
         case 'arfit2'
            [w,A_hat,Sigma_eps_hat,SBC] = arfit2(datChunk,pmin,pmax);
      end
      
%       if ~isempty(strfind(isnan(datChunk(:,1))',[1 1]))
%       keyboard;
%       end
         
      % find p if not fixed
      p = pmin;
      
      % MVAR_SPECTRAL expects the coefficient matrix to be 3-D, nxnxp
      A3D_hat = reshape(A_hat,n,n,p);
      mvar(i) = mvar_spectral(A3D_hat,Sigma_eps_hat,nfft,fs,'coherence',[],'granger',{l m});
      
      if size(datChunk,3) == 1
         ind = isnan(datChunk(:,1));
         datChunk(ind,:) = [];
      end
      
      [~,res] = arres(w,A_hat,datChunk,p+1);
      
      if nargout > 1
          % Breitung & Candelon significance test (2006)
          temp_res = reshape(res,size(res,1)*size(res,3),size(res,2));
          temp_datChunk = reshape(datChunk,size(datChunk,1)*size(datChunk,3),size(datChunk,2));
          bc(i) = freq_bc_test(omega_vec,[0 0]',A3D_hat,temp_datChunk,temp_res,l,m,p,alpha,fs);
      end
   end
end

if nargout > 2
    winStartT = (winStart - 1)*(1/fs);
end