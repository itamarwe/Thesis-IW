clear all;
close all;

C = 3e8;
Fc = 1e9;
Fs = 2^23;
N=512;
B = 0.4*Fs;
L = 4;

timeVec = (0:N-1)/Fs;

rTransmitter = [100,300];
v = [200,200];
rReceiverMat = LocateReceivers(1000,L,1);
[originalSignal, cleanReceivedSignal, sNorm2, sDotNorm2,imCorrssDot,delayedSignal, delayedSignalDot,b] = CreateRxSignal(rTransmitter, v, C, Fc, rReceiverMat,B,Fs,N);
%originalSignal = randn(1,512);
%originalSignal(200:300)=1;
%cleanReceivedSignal = [originalSignal;,circshift(originalSignal,[0,10]);circshift(originalSignal,[0,-5])];
%rReceiverMat = rReceiverMat(1:3,:);

%figure; plot(timeVec,abs(cleanReceivedSignal(1,:)),'b',timeVec, abs(cleanReceivedSignal(2,:)),'k');
%figure; plot(0:N-1,abs(fftshift(fft(cleanReceivedSignal(1,:)))),'b',0:N-1, abs(fftshift(fft(cleanReceivedSignal(2,:)))),'k');
%drawnow;
%profile on;
[DT,DF] = DT_DF_gridsearch(cleanReceivedSignal, rReceiverMat, rTransmitter,v,C,Fc,Fs);
%profile viewer;