% should add a case with jittered onset times of mean change, see Ding
% paper on artifacts

clear all; 

rng(10101); % fix random numbers

nTrials = 50;

flag = 4;

switch flag
   case 1
      % Constant underlying AR model with
      % mean change, different for both channels
      len = 1000;
      dat1 = sim_ar_model(1,len,nTrials);
      temp = [zeros(len,1) ,  [zeros(len/2,1) ; 4*ones(len/2,1)]];
      temp([len/5:len/5+200]) = 10;
      for i = 1:nTrials
         dat(i).v = dat1.data(i).v + temp;
      end
   case 2
      % Constant underlying AR model with
      % mean change, common to both channels
      len = 1000;
      dat1 = sim_ar_model(1,len,nTrials);
      temp = [zeros(500,1) ; 10*ones(500,1)];
      for i = 1:nTrials
         dat(i).v = dat1.data(i).v + repmat(temp,[1 2]);
      end
   case 3
      % Two different segments together from different AR models
      len = 500;
      dat1 = sim_ar_model(2,len,nTrials);
      dat2 = sim_ar_model(4,len,nTrials);
      
      for i = 1:nTrials
         dat(i).v = [dat1.data(i).v ; dat2.data(i).v];
      end
   case 4
      % Same AR model, different trial lengths. This only works with ARFIT2
      len = 100 + geornd(.25,nTrials,1)*20;
      dat1 = sim_ar_model(1,len,nTrials);
      for i = 1:nTrials
         dat(i).v = dat1.data(i).v;
         %dat(i).v(:,1) = dat(i).v(:,1)*3;
      end
   case 5
      % Same AR model, different trial lengths. This only works with ARFIT2
      % Add a mean change when the trial lengths start changing
      baselen = 100;
      len = baselen + geornd(.25,nTrials,1)*20;
      dat1 = sim_ar_model(2,len,nTrials);
      for i = 1:nTrials
         temp = [zeros(baselen,1) ; 10*ones(len(i)-baselen,1)];
         dat(i).v = dat1.data(i).v + repmat(temp,[1 2]);
      end
   otherwise
      error('Inappropriate flag');
end

l = 1;
m = 1;
winLen = 75;
stepLen = 15;
fit = 'arfit2';
pmin = 3;
pmax = 3;

fs = 1000;  % sampling frequency
%nfft = fs;  % length of fft
%demean = true;
%[mvar,bc,winStart] = sliding_mvar_spectral(dat,l,m,winLen,stepLen,fit,pmin,pmax,nfft,fs,demean);
[mvar,bc,winStartT] = sliding_mvar_spectral(dat,l,m,fs,winLen,stepLen,pmin,pmax);

f = mvar(1).f(1:fs/2);
t = winStart/1000;
% Power spectrum
S_y = [];
for i = 1:length(mvar)
   S_y = [S_y , abs(squeeze(mvar(i).S(1,1,1:fs/2)))];
end
S_z = [];
for i = 1:length(mvar)
   S_z = [S_z , abs(squeeze(mvar(i).S(2,2,1:fs/2)))];
end
% Coherence
C_yz = [];
for i = 1:length(mvar)
   C_yz = [C_yz , abs(squeeze(mvar(i).R(1,2,1:fs/2)))];
end
% Granger
G_yz = cat(2,mvar.G_yz);
G_yz = G_yz(1:fs/2,:);
G_zy = cat(2,mvar.G_zy);
G_zy = G_zy(1:fs/2,:);
p_yz = cat(2,bc.pval_yz);
p_zy = cat(2,bc.pval_zy);

figure;
subplot(511);
imagesc(t,f,log(S_y)); colorbar
title('Power y')
subplot(512);
imagesc(t,f,log(S_z)); colorbar
title('Power z')
subplot(513);
imagesc(t,f,C_yz); colorbar
title('Coherence yz')
subplot(514);
imagesc(t,f,G_yz); colorbar
title('Granger y<-z')
subplot(515);
imagesc(t,f,G_zy); colorbar
title('Granger z<-y')

figure;
subplot(511);
n = histc(len/1000,t);
bar(t,n,'histc'); 
axis([t(1) t(end) get(gca,'ylim')]);
colorbar
title('Trial lengths');
subplot(512);
imagesc(t,f,G_yz); colorbar
title('Granger y<-z')
subplot(513);
imagesc(t,f,p_yz<0.01); colormap gray; colorbar
title('Significance BC test')
subplot(514);
imagesc(t,f,G_zy); colorbar
title('Granger z<-y')
subplot(515);
imagesc(t,f,p_zy<0.01); colormap gray; colorbar
title('Significance BC test')


