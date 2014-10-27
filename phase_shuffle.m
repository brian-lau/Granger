% PHASE_SHUFFLE              Randomize phase of signal
% 
%     y = phase_shuffle(x);
%
%     INPUTS
%     x - input signal
%     OUTPUTS
%     y - phase-randomized signal

%     $ Copyright (C) 2014 Brian Lau http://www.subcortex.net/ $
%     Released under the BSD license. The license and most recent version
%     of the code can be found on GitHub:
%     https://github.com/brian-lau/Granger
%
%     REVISION HISTORY:
%     brian 11.14.07 written 

function y = phase_shuffle(x)

n = length(x);

if rem(n,2) == 0
   n2 = n/2;
else
   n2 = (n-1)/2;
end
fx = fft(x,2*n2);
magx = abs(fx);
phax = angle(fx);
phay = rand(n2-1,1) * 2*pi;
phay = [0; phay; phax(n2+1); -flipud(phay)];
% New Fourier transformed data with only the phase changed
magy = [magx(1:n2+1)' flipud(magx(2:n2))']';
fy = magy .* exp(phay.*i); 
% Transform back to time domain
y = real(ifft(fy,n));
