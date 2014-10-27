len = 1000;    % samples per signal
nsim = 100;    % # of simulations to run

l = 1;         % size of first subsystem
m = 1;         % size of second subsystem
n = l+m;       % total # of signals
pmin = 3;      % minimum AR order to estimate
pmax = 10;     % maximum AR order to estimate

fs = 1000;     % sampling frequency
nfft = fs;     % length of fft
alpha = 0.05;  % critical level for statistical tests

trueflag = logical(0); % use true generating coefficients, 0 forces a fit

model = 5;
omega = (pi)/2;% used for the simulations, in the models below, causality from 1<-2 will NOT
               % exist at this frequency, specify in NORMALIZED ANGULAR FREQUENCY

if model == 1
   %----- Basic system, Z granger causes Y, except at frequency OMEGA, page 369 Breitung & Candelon
   alpha = 0.3;
   A1 = [ 0.0   alpha;
          0.0   0.0 ];
   A2 = [ 0.0  -2*alpha*cos(omega);
          0.0  0.0 ];
   A3 = [ 0.0  alpha;
          0.0  0.0 ];
   Sigma_eps = [ 0.50  0.20;
                 0.20  0.50 ];
elseif model == 2
   %----- Simple low pass causality from Y<-Z
   A1 = [ 0.4   0.0;
          0.0   0.7 ];
   A2 = [ 0.35  0.2;
          0.0  -0.5 ];
   A3 = [ 0.0  0.1;
          0.0  0.0 ];
   Sigma_eps = [ 1.00  0.10;
                 0.10  0.50 ];
elseif model == 3
   %----- Basic system, Z granger causes Y, except at frequency OMEGA, page 372 Breitung & Candelon
   alpha = 0.3;
   A1 = [ 0.1   alpha;
          -1   0.1 ];
   A2 = [ 0.0  -2*alpha*cos(omega);
          0.0  -0.2 ];
   A3 = [ 0.0  alpha;
          0.0  0.3 ];
   Sigma_eps = [ 0.50  0.20;
                 0.20  0.50 ];
elseif model == 4
   % %----- Basic system, Z granger causes Y, except at frequency OMEGA, page 372 Breitung & Candelon
   % %----- also there is causality from Z<-Y of the lowpass form
   alpha = 0.3;
   A1 = [ 0.1   alpha;
          0.0   0.0 ];
   A2 = [ 0.0  -2*alpha*cos(omega);
          0.5  -0.0 ];
   A3 = [ 0.0  alpha;
          0.1  0.0 ];
   Sigma_eps = [ 1.00  0.10;
                 0.10  1.50 ];
elseif model == 5
   % %----- Basic system, Z granger causes Y, except at frequency OMEGA, page 372 Breitung & Candelon
   % %----- also there is causality from Z<-Y of the Gegenbaur form
   alpha = 0.3;
   A1 = [ 0.1   alpha;
          alpha   0.1 ];
   A2 = [ 0.0  -2*alpha*cos(omega);
          -2*alpha*cos(omega/2)  -0.0 ];
   A3 = [ 0.0  alpha;
          alpha  0.0 ];
   Sigma_eps = [ 0.50  0.20;
                 0.20  0.50 ];
else
   error('Model does not exist');
end

A = [ A1 A2 A3 ];
for i = 1:nsim

   v = arsim([0 0],A,Sigma_eps,len);

   if ~trueflag
      % Estimate Granger and BC stats using estimated model
      [w,A_hat,Sigma_eps_hat,SBC,FPE,th] = arfit(v,pmin,pmax);
      %[w,A_hat,Sigma_eps_hat,SBC,AIC] = arfit2(v,pmin,pmax);
      [siglev,res] = arres([0 0]',A_hat,v);
      ind = find(SBC==min(SBC));
      p = pmin:pmax;
      p = p(ind);
      A3D = reshape(A_hat,n,n,p);
   else
      % Estimate Granger and BC stats using true model
      p = length(A(:))/(n*n);
      [siglev,res] = arres([0 0]',A,v);
      A3D = reshape(A,n,n,p);
      Sigma_eps_hat = Sigma_eps;
   end

   mvar_point(i) = mvar_spectral(A3D,Sigma_eps_hat,nfft,fs,'granger',{l m});
   f = fftshift(mvar_point(1).f);
   f = mvar_point(1).f(1:nfft/2);

   omega_test = f;
   r = [0 0]'; % Test that magnitude is 0 for both COS and SIN components
   bc(i) = freq_bc_test(omega_test,r,A3D,v,res,l,m,p,alpha,fs);
end

if nsim > 1
   figure; 
   subplot(321); hold on
   f = mvar_point(1).f(1:nfft/2);
   temp = cat(2,mvar_point.G_yz);
   temp = temp(1:nfft/2,:);
   plot(f,temp,'-','Color',[.7 .7 .7]);
   plot(f,mean(temp,2),'k-','linewidth',3);
   ylim = get(gca,'ylim');
   axis([0 fs/2 0 ylim(2)]);
   title('Granger Y<-Z')
   
   subplot(322); hold on
   temp = cat(2,mvar_point.G_zy);
   temp = temp(1:nfft/2,:);
   plot(f,temp,'-','Color',[1 .75 .75]);
   plot(f,mean(temp,2),'r-','linewidth',3);
   ylim2 = get(gca,'ylim');
   if ylim2(2) > ylim(2)
      ylim = ylim2;
      subplot(321); axis([0 fs/2 0 ylim(2)]);
   end
   subplot(322)
   axis([0 fs/2 ylim]);
   title('Granger Z<-Y')
   
   subplot(323); hold on
   f = omega_test(:);
   temp = cat(2,bc.F_yz);
   plot(f,temp,'-','Color',[.7 .7 .7]);
   plot(f,mean(temp,2),'k-','linewidth',3);
   ylim = get(gca,'ylim');
   axis([0 fs/2 0 ylim(2)]);
   title('Breitung-Candelon Y<-Z')
   
   subplot(324); hold on
   temp = cat(2,bc.F_zy);
   plot(f,temp,'-','Color',[1 .75 .75]);
   plot(f,mean(temp,2),'r-','linewidth',3);
   ylim2 = get(gca,'ylim');
   if ylim2(2) > ylim(2)
      ylim = ylim2;
      subplot(323); axis([0 fs/2 0 ylim(2)]);
   end
   subplot(324);
   axis([0 fs/2 ylim]);
   title('Breitung-Candelon Z<-Y')
   
   subplot(325); hold on
   temp = cat(2,bc.pval_yz);
   temp = sum(temp<0.05,2)./nsim;
   plot(f,temp,'ko');
   plot([0 fs/2],[0.05 0.05],'b')
   axis([0 fs/2 0 1]);
   title('Breitung-Candelon Y<-Z reject at \alpha = 0.05')

   subplot(326); hold on
   temp = cat(2,bc.pval_zy);
   temp = sum(temp<0.05,2)./nsim;
   plot(f,temp,'ro');
   plot([0 fs/2],[0.05 0.05],'b')
   axis([0 fs/2 0 1]);
   title('Breitung-Candelon Z<-Y reject at \alpha = 0.05')
end