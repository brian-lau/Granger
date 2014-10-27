A1 = [ 0.4   0.7;
       0.0   0.7 ];
A2 = [ 0.35  -0.7;
       0.0  -0.5 ];
Sigma_eps = [ 1.00  0.10;
              0.10  1.50 ];
A = [ A1 A2 ];

v = arsim([0 0],A,Sigma_eps,200);

nfft = 256;
l = 1;
m = 1;
p = 2;

mvar = mvar_spectral(A,Sigma_eps,p,nfft,'coherence',[],'granger',{l m});

A1 = [ 0.5 0.5 0.0 0.0;
       -.0 0.5 0.0 0.0;
       0.0 0.0 0.5 0.5;
       0.0 0.0 .1 0.5];
A2 = [ 0.25 0.25 0.0 0.0;
       -.0 0.25 0.0 0.0;
       0.0 0.0 0.25 0.25;
       0.0 0.0 -.1 0.25];
% A1 = [ 0.4 0.31 0.3 0.3;
%      -.2 0.4 0.0 0.0;
%      0.0 0.0 0.4 0.0;
%      0.0 0.0 0.0 0.4];
% A2 = [ 0.3 0.5 0.1 0.6;
%      -.1 0.3 0.0 0.0;
%      0.0 0.0 0.3 0.0;
%      0.0 0.0 0.0 0.3];
Sigma_eps = [ 1.0 0.25 0.0 0.0;
              0.25 1.0 0.0 0.0;
              0.0 0.0 1.0 0.25;
              0.0 0.0 0.25 1.0];
A = [ A1 A2 ];

v = arsim([0 0 0 0],A,Sigma_eps,200);

nfft = 256;
l = 2;
m = 2;
p = 2;

mvar = mvar_spectral(A,Sigma_eps,p,nfft,'coherence',[],'granger',{l m});

%%%%%%%%%%%%%%
l = 2;
m = 2;
n = l+m;
p = 4;
filt1 = 0.5*(1-0.5).^(0:(p-1));
filt2 = 0.5*(1-0.5).^(0:(p-1));
filt3 = 0.5*(1-0.5).^(0:(p-1));
% filt2 = 0.5*(1-0.5).^(0:(p-1));
% filt3 = 0.75*(1-0.5).^(0:(p-1));

A = zeros(l+m,l+m,p);
A(1:l,1:l,:) = reshape(repmat(filt2,[l*l 1]),[l l p]);
A((l+1):n,(l+1):n,:) = reshape(repmat(filt3,[m*m 1]),[m m p]);
A(2,1,:) = 0;
A(4,3,:) = 0;

Sigma_eps = eye(n);%diag(1:n);
v = arsim(zeros(1,n),reshape(A,n,n*p),Sigma_eps,200);

% % A = [[0 0 ;0 0], [-.9 0; 0 0]]
% % C = eye(2);
% A1 = [ 0.4   1.2;   0.3   0.7 ];
% A2 = [ 0.35 -0.3;  -0.4  -0.5 ];
% % A1 = [ 0.4   0.7;
% %        0.0   0.7 ];
% % A2 = [ 0.35  -0.7;
% %        0.0  -0.5 ];
% Sigma_eps = [ 1.00  0.50;   0.50  1.50 ];
% A = [ A1 A2 ];
% 
% 
% 
% nfft = 256;
% l = 1;
% m = 1;
% %n = l + m;
% p = 2;
% 
% %[H,S,B] = mvar_spectral(A,Sigma_eps,p,nfft);
% mvar = mvar_spectral(A,Sigma_eps,p,nfft);
% [G_zy,G_yz] = fgranger(mvar.H,mvar.S,mvar.Sigma_eps,l,m);
% 
% % lind = 1:l;
% % mind = (l+1):(l+m);
% % 
% % % Partition the noise covariance
% % Sigma_eta = Sigma_eps(lind,lind);
% % Sigma_nu = Sigma_eps(mind,mind);
% % C = Sigma_eps(lind,mind);
% % 
% % B = [eye(n) , -A];
% % %A3D = reshape(B,l+m,l+m,p+1);
% % 
% % fB = fft(reshape(B,n,n,p+1),nfft,3);
% % 
% % H = zeros(size(fB));
% % for i = 1:nfft
% %    H(:,:,i) = inv(fB(:,:,i));
% % end
% % 
% % S = zeros(size(H));
% % for i = 1:nfft
% %    S(:,:,i) = H(:,:,i)*Sigma_eps*H(:,:,i)';
% % end
% % 
% % H_tilde = H;
% % H_tilde(lind,lind,:) = H(lind,lind,:) + H(lind,mind,:)*C'*inv(Sigma_eta);
% % H_tilde(mind,lind,:) = H(mind,lind,:) + H(mind,mind,:)*C'*inv(Sigma_eta);
% % 
% % G_zy = zeros(nfft,1);
% % for i = 1:nfft
% %    G_zy(i) = log( det(S(lind,lind,i)) / det(H_tilde(lind,lind,i)*Sigma_eta*H_tilde(lind,lind,i)') );
% % end
% % 
% % H_hat = H;
% % H_hat(lind,mind,:) = H(lind,mind,:) + H(lind,lind,:)*C*inv(Sigma_nu);
% % H_hat(mind,mind,:) = H(mind,mind,:) + H(mind,lind,:)*C*inv(Sigma_nu);
% % 
% % G_yz = zeros(nfft,1);
% % for i = 1:nfft
% %    G_yz(i) = log( det(S(mind,mind,i)) / det(H_hat(mind,mind,i)*Sigma_nu*H_hat(mind,mind,i)') );
% % end