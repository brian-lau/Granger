% Frequency vector for spectral analyses

function f = freq(fs,nfft)

deltaF = fs/nfft;
if rem(nfft,2), % ODD, don't include Nyquist.
   wmin = -(fs - deltaF)/2;
   wmax = (fs - deltaF)/2;   
else            % EVEN include Nyquist point in the negative freq.
   wmin = -fs/2;
   wmax = fs/2 - deltaF;
end
f = linspace(wmin, wmax, nfft);
if rem(nfft, 2) % ODD
   f((nfft+1)/2) = 0;
else
   f(nfft/2+1) = 0;
end

f = fftshift(f);