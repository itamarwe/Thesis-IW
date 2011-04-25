%TDOA FDOA Estimation test
close all;
clear all;

%Number of samples
N = 512;
%Pulse width
W = 30;
%Sampling Frequency
Fs = 512;

%Time delay and doppler
tau = 13/Fs;
m = floor(tau*Fs);
nu = -31; %[hz]

timeVec = ((0:N-1)/Fs)';

%Original Signal
x_n = zeros(N,1);
x_n((N/2-W/2):(N/2+W/2)) = 1;

%Received signal in receiver 1
y1_n = x_n+0.01*randn(size(x_n));

%Received signal in receiver 2
y2_n = circshift(x_n,-m).*exp(2*pi*j*nu*(timeVec+m))+0.01*randn(size(x_n));

%Plot the graphs
figure; plot(timeVec,x_n,'b',timeVec,y1_n,'r',timeVec, abs(y2_n),'k');
legend('x[n]','y1[n]','y2[n]');

TDOAVec = -20:20; %Samples
FDOAVec = -200:200; %[Hz]

[TDOAMat,FDOAMat] = meshgrid(TDOAVec,FDOAVec);

R = zeros(size(TDOAMat));
k=1;
for l=1:length(TDOAVec)
    ak = fft(y1_n);
    bk  = ak.*(exp(2*pi*1i/N*(0:N-1)*TDOAVec(l)).');
    y1_n_tshifted = ifft(bk);
    for k=1:length(FDOAVec)
        
        %y1_n_tshifted = circshift(y1_n,-TDOAMat(k,l));
        y1_n_ftshifted = exp(2*pi*1i*FDOAVec(k)*timeVec).*y1_n_tshifted;
        R(k,l) = 1\N*abs(y1_n_ftshifted'*y2_n);
    end
end


R = abs(R)/N;

figure;
mesh(TDOAMat,FDOAMat,abs(R));
xlabel('Time Shift [Samples]');
ylabel('Frequency Shift [Hz]');
