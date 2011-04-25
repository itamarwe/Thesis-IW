%profile on -timer real;
%1-gets a scenario from a file
%2-gets the parameters
%3-runs the estimations
%4-saves to the file

%****************************
%Get Parameters from file
%*******************************
[file dir] = uigetfile('experiment1.xls');
scenarios = xlsread([dir file]);

scenarioResults = zeros(size(scenarios,1),12); %[x,y,vx,vy,RMSXY, RMSVxVy, crbx,crby,crbvx,crbvy, CRBXY, CRBVxVy]

for scenario = 1:size(scenarios,1)
    C = 3e08; %Speed of Signal Propagation (Usually speed of light)[m/s]
    N = scenarios(scenario,1);
    L = scenarios(scenario,2); %Number of receivers
    Fs = scenarios(scenario,3); %Sampling Frequency [Hz]
%    TRIALS = scenarios(scenario,4); %Number of Monte Carlo Trials
    R = scenarios(scenario,4);% [m] Receivers Distance From the axis center


    %Transmitter Position, Velocity, Nominal Frequency, Bandwidth
    rTransmitter = [scenarios(scenario,5) scenarios(scenario,6)]; %Transmitter position [x,y] [m]
    v = [scenarios(scenario,7) scenarios(scenario,8)];%Transmitter velocity [vx, vy] [m/s]
    Fc = scenarios(scenario,9); % Carrier Frequency[Hz]
    B = scenarios(scenario,10);

    SNRdB = scenarios(scenario,11);
    receiverGeometry = scenarios(scenario,12);
    method = scenarios(scenario,13);
    estimations = scenarios(scenario,14);

    %*************************************
    % Run Experiment
    %*************************************

    %*** Receivers Positions ***
    %1 -Circular array with radius R[m] 2- Linear With apparatus 2R[m]
    rReceiverMat = LocateReceivers(R,L,receiverGeometry);

    estimationResults = zeros(estimations,4);

    for estimation = 1:estimations

        tic
        %gridResolution[x,y,vx,vy]
        %gridResolution  = [0.1 0.1 0.1 0.1];
        %gridPoints = floor(TRIALS^(0.25));

        %Calculate the longest delay for all of the points in the grid
        %maxGridSampleDelay = CalcMaxGridSampleDelay2 (rTransmitter(1),rTransmitter(2), gridResolution, gridPoints, rReceiverMat,C, Fs);
        maxGridSampleDelay = ceil(2*R*(Fs/C));
        %Create the original and received signal
        [originalSignal, cleanReceivedSignal, sNorm2, sDotNorm2 ,imCorrssDot,delayedSignal, delayedSignalDot,b] = ...
            CreateRxSignal(rTransmitter, v, C, Fc, rReceiverMat,B,Fs, N+maxGridSampleDelay);
        receivedSignal = addNoise(cleanReceivedSignal,SNRdB,N);

        disp(['Experiment No. ' num2str(estimation) '/' num2str(estimations)]);
        disp(['In Scenario ' num2str(scenario) '/' num2str(size(scenarios,1))]);

        %Perform Grid Search
        pScale =  20;
        vScale = 600;
        gridX = [linspace(-pScale/2, pScale/2,10)]+rTransmitter(1)+ (rand-.5)*10;
        gridY = [linspace(-pScale/2, pScale/2,10)]+rTransmitter(2)+(rand-.5)*10;
        gridVx = [linspace(-vScale/2, vScale/2,10)]+v(1)+(rand-.5)*10;
        gridVy = [linspace(-vScale/2, vScale/2,10)]+v(2)+(rand-.5)*10;

        while pScale > 0.005

            [x y vx vy] = GridSearch5(gridX, gridY, gridVx, gridVy, originalSignal, rReceiverMat, receivedSignal, Fc, C, Fs, N, method);
            %                [EmitPosEst_DPD,Xest_DPD,Yest_DPD]  =
            %                DPD_Known_Signals(time_vec,Xsearch_DPD,Ysearch_DPD,ReceiverPos,ReceiverVel,K,L,Q,fc,c,r_vec_t,CF1,sig);
            pScale=pScale*.6;
            vScale=vScale*.6;
            gridX = [0 linspace(-pScale/2, pScale/2,3)] + x;
            gridY = [0 linspace(-pScale/2, pScale/2,3)] + y;
            gridVx = [0 linspace(-vScale/2, vScale/2,3)]+vx;
            gridVy = [0 linspace(-vScale/2, vScale/2,3)]+vy;
        end

        % [x y vx vy] = GridSearch2(rTransmitter(1),rTransmitter(2),v(1),v(2),gridResolution,gridPoints, originalSignal, rReceiverMat, receivedSignal, Fc, C, Fs, N, method);
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
   % scenarioResults(scenario,7:10) = [sqrt(CRBMat(1,1)) sqrt(CRBMat(2,2)) sqrt(CRBMat(3,3)) sqrt(CRBMat(4,4))];
   % scenarioResults(scenario,11) = sqrt(scenarioResults(scenario,7)^2+scenarioResults(scenario,8)^2);
  %  scenarioResults(scenario,12) = sqrt(scenarioResults(scenario,9)^2+scenarioResults(scenario,10)^2);
    
    


    xlswrite([dir file '-results.xlsx'],[scenarios(1:scenario,:),scenarioResults(1:scenario,:)]);
end

%profile viewer;