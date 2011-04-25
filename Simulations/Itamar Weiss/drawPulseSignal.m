%Draw Pulse Signal

clear all;
close all;

C = 3e08; %Speed of Signal Propagation (Usually speed of light)[m/s]
N = 512;
L = 6; %Number of receivers
Fs = 2^23; %Sampling Frequency [Hz]
R = 1000;% [m] Receivers Distance From the axis center


%Transmitter Position, Velocity, Nominal Frequency, Bandwidth
rTransmitter = [100,300]; %Transmitter position [x,y] [m]
v = [200,200];%Transmitter velocity [vx, vy] [m/s]
Fc = 1e9; % Carrier Frequency[Hz]
B = 0.4*Fs;
receiverGeometry = 1; %1 - Circular 2- Linear

estimations = 10; %The number of estimations for each scenario


%Receivers Positions
rReceiverMat = LocateReceivers(R,L,receiverGeometry);

%The maximum sample delay because of the grid
maxGridSampleDelay = ceil(2*R*(Fs/C));

%Creating the received signal
[originalSignal, cleanReceivedSignal, sNorm2, sDotNorm2 ,imCorrssDot,delayedSignal, delayedSignalDot,b] = ...
    CreateRxSignal(rTransmitter, v, C, Fc, rReceiverMat,B,Fs, N+maxGridSampleDelay);

semilogy(linspace(-Fs/2,Fs/2,length(originalSignal)),abs(fftshift(fft(originalSignal))));
grid on;
title('Pulse Signal Power Spectrum','interpreter','latex','fontsize',16);
xtext='Frequency [Hz]';
xlabel(xtext,'interpreter','latex','fontsize',16);
ytext='Power [Arb]';
ylabel(ytext,'interpreter','latex','fontsize',16);
set(gca,'xtick', linspace(-Fs/2,Fs/2,5));
%set(gca,'xticklabel',linspace(-Fs/2,Fs/2,10))
figure;
plot((1:length(originalSignal))/Fs,abs(originalSignal));
grid on;
title('Pulse Signal','interpreter','latex','fontsize',16);
xtext='Time [Sec]';
xlabel(xtext,'interpreter','latex','fontsize',16);
ytext='Amplitude [Arb]';
ylabel(ytext,'interpreter','latex','fontsize',16);
%set(gca,'xtick', linspace(-Fs/2,Fs/2,5));
%set(gca,'xticklabel',linspace(-Fs/2,Fs/2,10))
