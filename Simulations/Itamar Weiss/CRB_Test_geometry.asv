
%In this script, we will calculate the CRB for an array of points
%in order to understand the relation between the performance and the
%geometry of the scenario.

%Script Structure:
%1 - Create the scenario (geometry, signals etc.)
%2 - Calculate numerical CRB for every point
%4 - Plots

close all;
clear all;

C = 3e08; %Speed of Signal Propagation (Usually speed of light)[m/s]
N = 512; %Number of smaples
L = 6; %Number of receivers
Fs = 2^15; %Sampling Frequency [Hz]
R = 1000;% [m] Receivers Distance From the axis center


%Transmitter Position, Velocity, Nominal Frequency, Bandwidth
%Transmitter examined positions
gridX = linspace(-4*R,4*R,100);
gridY = linspace(-4*R,4*R,100);
[gridX gridY] = meshgrid(gridX, gridY);

v = [200,200];%Transmitter velocity [vx, vy] [m/s]
Fc = 1e9; % Carrier Frequency[Hz]
B = 0.4*Fs;%Bandwidth[Hz]

SNRdB = 25; %The signal's SNR in all of the receivers
receiverGeometry = 1 %1-circular with radius R 2- Linear with length 2R

rReceiverMat = LocateReceivers(R,L,receiverGeometry); %The locations of the receivers

%%For a random array
%rReceiverMat = R*(randn(L,2)-0.5);

maxGridSampleDelay = ceil(8*R*(Fs/C)); %The maximum possible delay

dx = 1; %[m] the deltaX for the derivative
dy = 1;% [m] the delta y for the derivative
dvx = 1;%[m/s] the delta v_x for the derivative
dvy = 1; %[m/s] the delta v_y for the derivative

for i=1:size(gridX,1)
    i
    tic
    for j=1:size(gridY,2)

        rTransmitter = [gridX(i,j) gridY(i,j)];
%Create the original and received signal
[originalSignal, cleanReceivedSignal, sNorm2, sDotNorm2 ,imCorrssDot,delayedSignal, delayedSignalDot,b] = ...
    CreateRxSignalRand(rTransmitter, v, C, Fc, rReceiverMat,B,Fs, N+maxGridSampleDelay);

%d_m_l_d_x calculation

%m_l for x+dx
[originalSignal, ml_plus, sNorm2, sDotNorm2 ,imCorrssDot,delayedSignal, delayedSignalDot,b] = ...
    CreateRxSignal(rTransmitter+[dx 0], v, C, Fc, rReceiverMat,B,Fs, N+maxGridSampleDelay);
%m_l for x-dx
[originalSignal, ml_minus, sNorm2, sDotNorm2 ,imCorrssDot,delayedSignal, delayedSignalDot,b] = ...
    CreateRxSignal(rTransmitter+[-dx 0], v, C, Fc, rReceiverMat,B,Fs, N+maxGridSampleDelay);

%dMldX = (ml_plus(:,end-N+1:end)-ml_minus(:,end-N+1:end))/(2*dx);
dMldX = (ml_plus(:,1:N)-ml_minus(:,1:N))/(2*dx);

%d_m_l_d_y calculation

%m_l for y+dy
[originalSignal, ml_plus, sNorm2, sDotNorm2 ,imCorrssDot,delayedSignal, delayedSignalDot,b] = ...
    CreateRxSignal(rTransmitter+[0 dy], v, C, Fc, rReceiverMat,B,Fs, N+maxGridSampleDelay);
%m_l for x-dx
[originalSignal, ml_minus, sNorm2, sDotNorm2 ,imCorrssDot,delayedSignal, delayedSignalDot,b] = ...
    CreateRxSignal(rTransmitter+[0 -dy], v, C, Fc, rReceiverMat,B,Fs, N+maxGridSampleDelay);

%dMldY = (ml_plus(:,end-N+1:end)-ml_minus(:,end-N+1:end))/(2*dy);
dMldY = (ml_plus(:,1:N)-ml_minus(:,1:N))/(2*dy);

%d_m_l_d_vx calculation

%m_l for vx+dvx
[originalSignal, ml_plus, sNorm2, sDotNorm2 ,imCorrssDot,delayedSignal, delayedSignalDot,b] = ...
    CreateRxSignal(rTransmitter, v+[dvx 0], C, Fc, rReceiverMat,B,Fs, N+maxGridSampleDelay);
%m_l for vx-dvx
[originalSignal, ml_minus, sNorm2, sDotNorm2 ,imCorrssDot,delayedSignal, delayedSignalDot,b] = ...
    CreateRxSignal(rTransmitter, v+[-dvx 0], C, Fc, rReceiverMat,B,Fs, N+maxGridSampleDelay);

%dMldVx = (ml_plus(:,end-N+1:end)-ml_minus(:,end-N+1:end))/(2*dvx);
dMldVx = (ml_plus(:,1:N)-ml_minus(:,1:N))/(2*dvx);

%d_m_l_d_vy calculation

%m_l for vy+dvy
[originalSignal, ml_plus, sNorm2, sDotNorm2 ,imCorrssDot,delayedSignal, delayedSignalDot,b] = ...
    CreateRxSignal(rTransmitter, v+[0 dvy], C, Fc, rReceiverMat,B,Fs, N+maxGridSampleDelay);
%m_l for vy-dvy
[originalSignal, ml_minus, sNorm2, sDotNorm2 ,imCorrssDot,delayedSignal, delayedSignalDot,b] = ...
    CreateRxSignal(rTransmitter, v+[0 -dvy], C, Fc, rReceiverMat,B,Fs, N+maxGridSampleDelay);

%dMldVy = (ml_plus(:,end-N+1:end)-ml_minus(:,end-N+1:end))/(2*dvy);
dMldVy = (ml_plus(:,1:N)-ml_minus(:,1:N))/(2*dvy);

dMldX = dMldX';
dMldY = dMldY';
dMldVx = dMldVx';
dMldVy = dMldVy';

    
    sigmal = 10^(-SNRdB/20)*ones(1,L)/sqrt(N);

    FIM = zeros(4); %(1-x,2-y,3-vx,4-vy];
    %Calculation of the Fischer Information Matrix


    for m=1:4
        switch m
            case {1}
                Mlm=dMldX;
            case{2}
                Mlm = dMldY;
            case {3}
                Mlm = dMldVx;
            case {4}
                Mlm = dMldVy;
        end %switch m

        for n=1:4
            switch n
                case {1}
                    Mln=dMldX;
                case{2}
                    Mln = dMldY;
                case {3}
                    Mln = dMldVx;
                case {4}
                    Mln = dMldVy;
            end %switch n

            for l=1:L
                FIM(m,n)=FIM(m,n)+(2/(sigmal(l)^2))*real(Mlm(:,l)'*Mln(:,l));
            end %l
        end %n
    end%m

    CRBMat = inv(FIM);
    CRBX(i,j) = CRBMat(1,1);
    CRBY(i,j) = CRBMat(2,2);
    CRBVx(i,j) = CRBMat(3,3);
    CRBVy(i,j) = CRBMat(4,4);
    end%j
    toc
end%i

plot(rReceiverMat(:,1),rReceiverMat(:,2),'kp');
hold on;
contour(gridX, gridY, 0.5*log10(CRBX+CRBY));
title 'Log10 RMS Position Error';
colorbar;
xtext = 'X[m]';
ytext = 'Y[m]';
legend('Receivers');

figure;
plot(rReceiverMat(:,1),rReceiverMat(:,2),'kp');
hold on;
contour(gridX, gridY, 0.5*log10(CRBVx+CRBVy));
title 'RMS Velocity Error';
colorbar;
xtext = 'X[m]';
ytext = 'Y[m]';
legend('Receivers');