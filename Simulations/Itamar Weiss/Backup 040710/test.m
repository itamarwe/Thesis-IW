clear all;
close all;
N = 128;
W = 11;
k = 0:N-1;
Fs = 1;
B = 0.8*Fs;

%Fourier transform of a 1/W high step with width W
ak = sin(pi*W*k/N)./(W*sin(pi*k/N));
ak(1) = 1;

%delayed coeffients dx[n] = x[n-m];
m = floor(N/2)+0.2;
bk = ak.*exp(-2*pi*i/N*k*m);

%Reconstruction

tic
[k,n] = meshgrid(0:N-1,0:N-1);
F = exp(2*pi*i/N*k.*n);
xn = 1/N*ak*F;
dxn = 1/N*bk*F;
toc
tic
xn = ifft([ak 0]);
dxn = ifft([dxn 0]);
toc

%Band Limited Step
%lets say we want to keep K coefficients from each side

K= N/2-10;
ck = bk;
ck(K+2:N-K)=0;
%create delay for ck
bxn = 1/N*ck*F;

figure;stem(0:N-1,ak);title('a[k] coefficients');
figure; stem(0:N-1,xn);title('reconstructed x[n]');
hold on; stem(0:N-1,dxn,'r');stem(0:N-1,bxn,'g');

%making sure that using cosins gives the same result

dk = 2*bk;
dk(ceil(N/2):end)=0;
dk(1) = bk(1);
cosxn = 1/N*dk*cos(2*pi/N*k.*n);

figure;stem(0:N-1,cosxn);

figure; stem(real(fft(xn))+0.01); hold on; stem(real(ak),'r');
figure; stem(imag(xn)+0.01); hold on; stem(imag(ifft(ak)),'r');