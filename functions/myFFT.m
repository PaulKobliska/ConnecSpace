function [LFP_FFT] = myFFT(LFP,Fs);

% this function calculate the fast fourier transform of the LFP and plot it
% for the entore frequency range (i.e. 512 for a sampling frequency of 1024
% Hz, or the Nequist frequency)
% Fs = sampling freq (ex: 1024)
% LFP= amplitude LFP signal (1 column)


T = 1/Fs;         
L = length(LFP);
t = (0:L-1)*T; % time vector
NFFT = 2^nextpow2(L); % transform data points on power2 bases
FFT = fft(LFP,NFFT)/L;
f = Fs/2*linspace(0,1,NFFT/2+1);
f=f';
LFP_FFT=zeros(length(f),2);
LFP_FFT(:,1)=f;
LFP_FFT(:,2)=2*abs(FFT(1:NFFT/2+1));
