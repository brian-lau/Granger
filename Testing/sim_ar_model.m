% Simulate some simple systems for testing functions

function dat = sim_ar_model(model,len,nsim)

% model = 1;     % see below for the different model types
% len = 1000;    % samples per signal
% nsim = 10;     % # of simulations to run

if length(len) > 1
   if length(len) ~= nsim
      error('length of LEN must equal NSIM');
   end
else
   len = repmat(len,nsim,1);
end

omega = (pi)/2;% used for the simulations, in the models below, causality from 1<-2 will NOT
               % exist at this frequency, specify in NORMALIZED ANGULAR FREQUENCY

if model == 1
   % Basic system, causality from Y<-Z, except at frequency OMEGA, page 369 Breitung & Candelon
   l = 1;
   m = 1;
   beta = 0.3;
   A1 = [ 0.0   beta;
          0.0   0.0 ];
   A2 = [ 0.0  -2*beta*cos(omega);
          0.0  0.0 ];
   A3 = [ 0.0  beta;
          0.0  0.0 ];
   Sigma_eps = [ 0.50  0.20;
                 0.20  0.50 ];
elseif model == 2
   % Simple low pass causality from Y<-Z
   l = 1;
   m = 1;
   A1 = [ 0.4   0.0;
          0.0   0.7 ];
   A2 = [ 0.35  0.2;
          0.0  -0.5 ];
   A3 = [ 0.0  0.1;
          0.0  0.0 ];
   Sigma_eps = [ 1.00  0.10;
                 0.10  0.50 ];
elseif model == 3
   % Basic system, causality from Y<-Z, except at frequency OMEGA, page 372 Breitung & Candelon
   l = 1;
   m = 1;
   beta = 0.3;
   A1 = [ 0.1   beta;
          -1   0.1 ];
   A2 = [ 0.0  -2*beta*cos(omega);
          0.0  -0.2 ];
   A3 = [ 0.0  beta;
          0.0  0.3 ];
   Sigma_eps = [ 0.50  0.20;
                 0.20  0.50 ];
elseif model == 4
   % Basic system, causality from Y<-Z, except at frequency OMEGA, page 372 Breitung & Candelon
   % also there is causality from Z<-Y of the lowpass form
   l = 1;
   m = 1;
   beta = 0.3;
   A1 = [ 0.1   beta;
          0.0   0.0 ];
   A2 = [ 0.0  -2*beta*cos(omega);
          0.5  -0.0 ];
   A3 = [ 0.0  beta;
          0.1  0.0 ];
   Sigma_eps = [ 1.00  0.10;
                 0.10  1.50 ];
elseif model == 5
   % Basic system, causality from Y<-Z, except at frequency OMEGA, page 372 Breitung & Candelon
   % also there is causality from Z<-Y, except at frequency OMEGA/2
   l = 1;
   m = 1;
   beta = 0.3;
   A1 = [ 0.1   beta;
          beta   0.1 ];
   A2 = [ 0.0  -2*beta*cos(omega);
          -2*beta*cos(omega/2)  -0.0 ];
   A3 = [ 0.0  beta;
          beta  0.0 ];
   Sigma_eps = [ 0.50  0.20;
                 0.20  0.50 ];
else
   error('Model does not exist');
end

A = [ A1 A2 A3 ];

dat.model = model;
dat.A = A;
dat.l = l;
dat.m = m;
dat.n = l+m;
dat.p = 3;
dat.Sigma_eps = Sigma_eps;

for i = 1:nsim
   % Simulate AR process
   %v(:,:,i) = arsim([0 0],A,Sigma_eps,len);
   dat.data(i).v = arsim([0 0],A,Sigma_eps,len(i));
end

%dat.data = v;
