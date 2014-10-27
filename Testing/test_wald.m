A1 = [ 0.4   0.7;
       0.0   0.7 ];
A2 = [ 0.35  -0.7;
       0.0  -0.5 ];
Sigma_eps = [ 1.00  0.10;
              0.10  1.50 ];
A = [ A1 A2 ];

v = arsim([0 0],A,Sigma_eps,200);

l = 1;      % size of first set
m = 1;      % size of second set
n = l+m;    % total size of system
pmin = 1;
pmax = 20;

[w,A_hat,Sigma_eps_hat,SBC,FPE,th] = arfit(v,pmin,pmax);
[siglev,res] = arres([0 0]',A_hat,v);
A3D = reshape(A_hat,n,n,find(SBC==min(SBC)));
% mvar_point = mvar_spectral(A3D,Sigma_eps_hat,nfft,fs,'coherence',[],'granger',{l m});

p = 2;
X = v;
%X = [(11:19)',(21:29)'];
T = length(X);

Xdot_temp = zeros(T,n*p);
for i = 1:n
   ind = (i-1)*2 + (1:p);
   temp = convmtx(X(:,i),p);
   Xdot_temp(:,ind) = [zeros(1,p) ; temp(1:end-p,:)];
end

Lambda_temp = zeros(n*p,n);
for i = 1:n
   temp = A3D(i,:,:);
   Lambda_temp(:,i) = reshape(squeeze(temp)',n*p,1);
end

Xdot = zeros(T,n*p);
for i = 1:n
   ind1 = (i-1)*2 + (1:p);
   ind2 = i:p:(n*p);
   Xdot(:,ind1) = Xdot_temp(:,ind2);
end

Lambda = zeros(n*p,n);
for i = 1:n
   temp = A3D(i,:,:);
   Lambda(:,i) = temp(:);
end

%res_temp = X - (Lambda_temp'*Xdot_temp')';
res = X - (Lambda'*Xdot')';

Sigma_hat = (res'*res)/(T-n*p);

% Sigma_hat = zeros(n,n);
% for i = 1:T
%    Sigma_hat = Sigma_hat + res(i,:)'*res(i,:);
% end
% Sigma_hat = Sigma_hat / (T-n*p);

var_vecLambda = kron(Sigma_hat,inv(Xdot'*Xdot));

H = eye(8);
h = zeros(8,1);
wald = (H*Lambda(:) - h)'*inv(H*var_vecLambda*H')*(H*Lambda(:) - h)
1-chi2cdf(wald,rank(H))

pval = mvar_tgranger_test(X,A3D,H,h)

H = [0 0 0 0 1 0 0 0 ; 
     0 0 0 0 0 0 1 0];
h = [0 0]';
wald = (H*Lambda(:) - h)'*inv(H*var_vecLambda*H')*(H*Lambda(:) - h)
1-chi2cdf(wald,rank(H))

R{1} = zeros(n,n,p);
R{1}(2,1,1) = 1;
r{1} = 0;
R{2} = zeros(n,n,p);
R{2}(2,1,2) = 1;
r{2} = 0;
pval = mvar_tgranger_test(X,A3D,R,r)

H = [0 0 0 0 0 0 1 0];
h = [0]';
wald = (H*Lambda(:) - h)'*inv(H*var_vecLambda*H')*(H*Lambda(:) - h)
1-chi2cdf(wald,rank(H))

clear R r;
R{1} = zeros(n,n,p);
R{1}(2,1,2) = 1;
r{1} = 0;
pval = mvar_tgranger_test(X,A3D,R,r)

H = [0 1 0 0 0 -1 0 0];
h = [0]';
wald = (H*Lambda(:) - h)'*inv(H*var_vecLambda*H')*(H*Lambda(:) - h)
1-chi2cdf(wald,rank(H))

clear R r;
R{1} = zeros(n,n,p);
R{1}(1,2,1) = 1;
R{1}(2,2,1) = -1;
r{1} = 0;
pval = mvar_tgranger_test(X,A3D,R,r)
