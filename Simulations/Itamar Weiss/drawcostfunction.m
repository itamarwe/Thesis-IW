%This script creates a 2-d graph of the Cost function.

clear all;
close all;

C = 3e08; %Speed of Signal Propagation (Usually speed of light)[m/s]
N = 512; %The number of samples
L = 6; %Number of receivers
Fs = 2^18; %Sampling Frequency [Hz]
R = 1000;% [m] Receivers Distance From the axis center


%Transmitter Position, Velocity, Nominal Frequency, Bandwidth
rTransmitter = [100,300]; %Transmitter position [x,y] [m]
v = [200,200];%Transmitter velocity [vx, vy] [m/s]

xMin = 0;
xMax = 200;
yMin = 200;
yMax = 400;


Fc =1e9; % Carrier Frequency[Hz]
B = 0.4*Fs;

SNRdB = inf;
receiverGeometry = 1;
method = 2;

rTransmitterOriginal = rTransmitter;
%*************************************
% Run Experiment
%*************************************

%*** Receivers Positions ***
rReceiverMat = LocateReceivers(R,L,receiverGeometry);

CRB_toa =0;
CRB_foa = 0;
DT = [];
DF= [];

%Calculates the CRB necessary for determining the weights in the
%conventional method
if (method==5)
    maxGridSampleDelay = ceil(2*R*(Fs/C));
    [originalSignal, cleanReceivedSignal, sNorm2, sDotNorm2 ,imCorrssDot,delayedSignal, delayedSignalDot,b] = ...
        CreateRxSignal(rTransmitter, v, C, Fc, rReceiverMat,B,Fs, N+maxGridSampleDelay);
    [CRB_toa CRB_foa] = CRB_toa_foa(Fc,C, delayedSignal, delayedSignalDot, N, Fs, b, rTransmitter,v, rReceiverMat,SNRdB);
    CRB_toa = max(CRB_toa);
    CRB_foa = max(CRB_foa);
end


%Calculate the longest delay for all of the points in the grid
%maxGridSampleDelay = CalcMaxGridSampleDelay2 (rTransmitter(1),rTransmitter(2), gridResolution, gridPoints, rReceiverMat,C, Fs);
%%maxGridSampleDelay = ceil(4*R*(Fs/C));
maxGridSampleDelay=0;
%Create the original and received signal
[originalSignal, cleanReceivedSignal, sNorm2, sDotNorm2 ,imCorrssDot,delayedSignal, delayedSignalDot,b] = ...
    CreateRxSignal(rTransmitter, v, C, Fc, rReceiverMat,B,Fs, N+2*maxGridSampleDelay);
receivedSignal = addNoise(cleanReceivedSignal,SNRdB,N);

if (method==5)
    %Conventional FDOA-TDOA methods
    %Find the DT and DF matrices
    [DT,DF] = DT_DF_gridsearch(receivedSignal, rReceiverMat, rTransmitter,v,C,Fc,Fs);
end

Ntilde = size(receivedSignal,2);

%A vector used for time translations in the frequency domain
k=0:Ntilde-1;

if (method==5)
    %The weights for the cost function
    w_toa = 1/CRB_toa;
    w_foa = 1/CRB_foa;
end

%Perform Grid Search

gridX = [linspace(xMin, xMax,100)];
%gridX = 100;
gridY = [linspace(yMin, yMax,100)];
%gridY = 300;
gVx = v(1);
%gVx = linspace(180,220,100);
gVy = v(2);
%gVy = linspace(180,220,100);

%costMat = zeros(length(gridX),length(gridY));
costMat = zeros(length(gVx),length(gVy));
h = waitbar(0,'Please Wait');
for gVxindex = 1:length(gVx)
            %waitbar(gVxindex/length(gVx),h);
    for gVyindex = 1:length(gVy)
        for gXindex = 1:length(gridX)
            waitbar(gXindex/length(gridX),h);
            for gYindex = 1:length(gridY)

                %The examined position of the transmitter
                rTransmitter = [gridX(gXindex) gridY(gYindex)];
                %The examined velocity of the transmitter
                v = [gVx(gVxindex) gVy(gVyindex)];
                %A matrix used for the calculation of the distances
                rTransmitterMat = ones(L,2)*[rTransmitter(1) 0;0 rTransmitter(2)];
                %A matrix with the difference vectors
                rDiffMat = rTransmitterMat-rReceiverMat;
                %A vector holding the distances between each receiver and
                %the transmitter
                rDistances = sqrt(rDiffMat(:,1).^2+rDiffMat(:,2).^2);
                %A calculation of the delay and maximum delay of the signal
                [maxDelay delay] = CalcMicroSecRxSignalDelay(rTransmitter, rReceiverMat,C);

                if (method==1)||(method==2)
                    %The doppler shift of the signals
                    miu = zeros(L,1);
                    dopplerFixedSignal = zeros(L,size(receivedSignal,2));
                    V = zeros(N,L);
                    for l=1:L
                        %miu calculation
                        %miu(l) = -1/C*v*rDiffMat(l,:)'/rDistances(l);
                        hzDopplerShift = -(Fc/C)*v*(rDiffMat(l,:)'/rDistances(l));
                        %frequency calculation
                        %hzDopplerShift = Fc*miu(l);
                        %Al Matrix
                        %hermitianA(:,:,l) = diag(exp(2*pi*i*hzDopplerShift*(1:N)/Fs))';
                        %dopplerFixedSignal(l,:) = (hermitianA(:,:,l)*receivedSignal(l,:).').';
                        
                            %%vHermitianA = exp(-2*pi*1i*hzDopplerShift*(1:size(receivedSignal,2))/Fs);
                            %%dopplerFixedSignal(l,:) =(vHermitianA.*receivedSignal(l,:)).';
                        
                        %dopplerFixedSignal(l,:) = time_freq_shift(receivedSignal(l,:),Fs,0,-hzDopplerShift);
                        
                           %% m = delay(l)*(Fs*1e-6); %04/07/10
                          %%  vk = fft(dopplerFixedSignal(l,:));
                        % vk = vk.*exp(2*pi*i/Ntilde*k*sampleDelay(l));% 04/07/10
                       %% vk = vk.*exp(2*pi*1i/Ntilde*k*m); % 04/07/10
                       %% vl = ifft(vk);
                       vl= time_freq_shift(receivedSignal(l,:),Fs,-delay(l)*1e-6,-hzDopplerShift);
                        V(:,l) = vl(1:N);
                        %        V(:,l) = dopplerFixedSignal(l,(1+sampleDelay(l)):(N+sampleDelay(l)))';
                        % V(:,l) = dopplerFixedSignal(l,(1+sampleMaxDelay-sampleDelay(l)):(N+sampleMaxDelay-sampleDelay(l)))';
                    end
                end
                %Evaluate Qtilda
                if (method==1)||(method==3)
                    %DPD Unknown Signals
                    Qtilda = V'*V;
                    %Evaluate LambdaMax(Qtilda)
                    cost = max(eig(Qtilda));
                elseif (method==2)||(method==4)
                    %DPD Known Signals
                    cost = originalSignal(1:N)*V*(V'*originalSignal(1:N)');
                elseif (method==5)
                    %Conventional method
                    DT0 = zeros(L);
                    DF0 = zeros(L);
                    for k2=1:L
                        for l2 = (k2+1):L
                            %DT0 calculation
                            DT0(k2,l2) = norm((rReceiverMat(k2,:)-rTransmitter))/C-norm((rReceiverMat(l2,:)-rTransmitter))/C;
                            %DF0 calculation
                            DF0(k2,l2) = -Fc/C*((v*(rTransmitter-rReceiverMat(k2,:))')/norm((rTransmitter-rReceiverMat(k2,:))))+Fc/C*((v*(rTransmitter-rReceiverMat(l2,:))')/norm((rTransmitter-rReceiverMat(l2,:))));
                        end
                    end
                    %The cost function of the conventional method
                    %The minus sign is put here so that maximum of the cost
                    %function (with the minus) will be in the minimum of
                    %the cost function (without the minus)
                    cost = -(w_toa*sum(sum((DT-DT0).^2))+w_foa*sum(sum((DF-DF0).^2)));
                end

                costMat(gXindex,gYindex) = cost;
                    %costMat(gVxindex,gVyindex) = cost;
            end %gYindex
        end%gXindex
    end %gVyindex
end%gVxindex
        close(h);

        figure;
        plot(rTransmitterOriginal(1),rTransmitterOriginal(2),'kp');
        %plot(200,200,'kp');
        hold on;
        contour(gridX,gridY,costMat);
        %contour(gVx,gVy,costMat);
        legend('Transmitter');
        xlabel('x[m]');
        ylabel('y[m]');
