%In this experiment we will examine the performance Vs. Number of Receivers
clear all;
close all;

C = 3e08; %Speed of Signal Propagation (Usually speed of light)[m/s]
N = 512;
L = 6; %Number of receivers
Fs = 2^20; %Sampling Frequency [Hz]
R = 1000;% [m] Receivers Distance From the axis center


%Transmitter Position, Velocity, Nominal Frequency, Bandwidth
rTransmitter = [100,1000]; %Transmitter position [x,y] [m]
v = [200,200];%Transmitter velocity [vx, vy] [m/s]
Fc = 1e9; % Carrier Frequency[Hz]
B = 0.4*Fs;
receiverGeometry = 2; %1 - Circular 2- Linear

estimations = 10; %The number of estimations for each scenario


%Receivers Positions
rReceiverMat = LocateReceivers(R,L,receiverGeometry);

%The maximum sample delay because of the grid
maxGridSampleDelay = ceil(8*R*(Fs/C));

%Creating the received signal
[originalSignal, cleanReceivedSignal, sNorm2, sDotNorm2 ,imCorrssDot,delayedSignal, delayedSignalDot,b] = ...
    CreateRxSignalRand(rTransmitter, v, C, Fc, rReceiverMat,B,Fs, N+maxGridSampleDelay);

SNR=25:-5:-15;

scenarioResults = zeros(3*length(SNR),12); %[x,y,vx,vy,RMSXY, RMSVxVy, crbx,crby,crbvx,crbvy, CRBXY, CRBVxVy]
scenario = 1;
for method = [1,5]
    for SNRdB = SNR
        
        CRB_toa =0;
        CRB_foa = 0;
        DT = [];
        DF= [];

        %Calculates the CRB necessary for determining the weights in the
        %conventional method
        if (method==5)
            [CRB_toa CRB_foa] = CRB_toa_foa(Fc,C, delayedSignal, delayedSignalDot, N, Fs, b, rTransmitter,v, rReceiverMat,SNRdB);
            CRB_toa = max(CRB_toa);
            CRB_foa = max(CRB_foa);
        end

        estimationResults = zeros(estimations,4);

        for estimation = 1:estimations

            tic
            %Create the noised signal
            receivedSignal = addNoise(cleanReceivedSignal,SNRdB,N);

            disp(['Experiment No. ' num2str(estimation) '/' num2str(estimations)]);
            disp(['SNRdB: ' num2str(SNRdB) ' Method: ' num2str(method)]);

            if (method==5)
                %Conventional FDOA-TDOA methods
                %Find the DT and DF matrices
                [DT,DF] = DT_DF_gridsearch(receivedSignal, rReceiverMat, rTransmitter,v,C,Fc,Fs);
            end

            %Perform Grid Search
            pScale =  10;
            vScale = 1000;
            gridX = linspace(-pScale/2, pScale/2,10)+rTransmitter(1)+ (rand-.5)*pScale/5;
            gridY = linspace(-pScale/2, pScale/2,10)+rTransmitter(2)+(rand-.5)*pScale/5;
            gridVx = linspace(-vScale/2, vScale/2,10)+v(1)+(rand-.5)*vScale/5;
            gridVy = linspace(-vScale/2, vScale/2,10)+v(2)+(rand-.5)*vScale/5;

            while pScale > 0.01

                [x y vx vy] = GridSearch5(gridX, gridY, gridVx, gridVy, originalSignal, rReceiverMat, receivedSignal, Fc, C, Fs, N, method,DT,DF,CRB_toa, CRB_foa);
                pScale=pScale*.6;
                vScale=vScale*.6;
                gridX = [linspace(-pScale/2, pScale/2,3)] + x;
                gridY = [linspace(-pScale/2, pScale/2,3)] + y;
                gridVx = [linspace(-vScale/2, vScale/2,3)]+vx;
                gridVy = [linspace(-vScale/2, vScale/2,3)]+vy;
            end
            disp('DPD estimation:');
            disp(['X:' num2str(x) ' Y:' num2str(y) ' Vx:' num2str(vx) ' Vy:' num2str(vy)]);
            disp(' ');

            estimationResults(estimation,:) = [x y vx vy];
            estimationResults(estimation,:) = estimationResults(estimation,:)-[rTransmitter v];
            toc
        end

        %***********************************
        %Cramer-Rao Lower Bound Calculation
        %CRBMat = CRB2(Fc,C, delayedSignal, delayedSignalDot, N, Fs, b, rTransmitter,v, rReceiverMat,SNRdB);
        %***********************************

        scenarioResults(scenario,1:4) = sqrt(mean(estimationResults.^2));
        scenarioResults(scenario,5) = sqrt(scenarioResults(scenario,1)^2+scenarioResults(scenario,2)^2);
        scenarioResults(scenario,6) = sqrt(scenarioResults(scenario,3)^2+scenarioResults(scenario,4)^2);
%        scenarioResults(scenario,7:10) = [sqrt(CRBMat(1,1)) sqrt(CRBMat(2,2)) sqrt(CRBMat(3,3)) sqrt(CRBMat(4,4))];
  %      scenarioResults(scenario,11) = sqrt(scenarioResults(scenario,7)^2+scenarioResults(scenario,8)^2);
   %     scenarioResults(scenario,12) = sqrt(scenarioResults(scenario,9)^2+scenarioResults(scenario,10)^2);
        xlswrite('experiment5-results4.xlsx',[N,L,Fs,R,rTransmitter(1),rTransmitter(2),v(1),v(2),Fc,B,SNRdB,receiverGeometry,method,estimations,...
            scenarioResults(scenario,:)],['A' num2str(scenario) ':Z' num2str(scenario)]);
        scenario = scenario+1;
    end
end





%profile viewer;