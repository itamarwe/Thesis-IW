function [x,y,vx,vy] = GridSearch2(cX,cY,cVx,cVy,gridResolution, gridPoints, originalSignal, rReceiverMat, receivedSignal, Fc, C, Fs, N, method)
%function [x,y,Vx,Vy] = GridSearch2(cX,cY,cVx,cVy,gridReolution,gridPoints, originalSignal, rReceiverMat, receivedSignal, Fc, C, Fs, N, method)
%cX,cY,cVx,cVy - the center of the grid
%gridResolution - 1x4 vector with the resolution in every parameter
%gridPoint - number of grid points for every parameter
% method: 1 - UnknownSignals, 2 - KnownSignals, 3 - KnownNarrowBandSignals
%**********************************************************8
% GRID SEARCH *********************************8
% *************************************************8
L = size(rReceiverMat,1);
maxCost = - inf;
waitbarHandle = waitbar(0,'Processing');
trial=0;
TRIALS = gridPoints^4;
gridX = ((-(gridPoints-1)/2):1:((gridPoints-1)/2))*gridResolution(1)+cX;
gridY = ((-(gridPoints-1)/2):1:((gridPoints-1)/2))*gridResolution(2)+cY;
gridVx = ((-(gridPoints-1)/2):1:((gridPoints-1)/2))*gridResolution(3)+cVx;
gridVy = ((-(gridPoints-1)/2):1:((gridPoints-1)/2))*gridResolution(4)+cVy;

for gX = gridX
    for gY = gridY
        for gVx = gridVx
            for gVy = gridVy
                trial = trial+1;
                if mod(trial,10)==0
                    waitbar(trial/TRIALS,waitbarHandle);
                end

                rTransmitter = [gX gY];
                v = [gVx gVy];

                rTransmitterMat = ones(L,2)*[rTransmitter(1) 0;0 rTransmitter(2)];
                rDiffMat = rTransmitterMat-rReceiverMat;
                rDistances = sqrt(rDiffMat(:,1).^2+rDiffMat(:,2).^2);
                microsecDelay = rDistances/(C/1e6);
                sampleDelay = floor(microsecDelay*Fs/1e6);
                sampleMaxDelay = max(sampleDelay);
                %Evaluate hermitian A
                %    hermitianA = zeros(N,N,L);
                miu = zeros(L,1);
                dopplerFixedSignal = zeros(L,size(receivedSignal,2));
                for l=1:L
                    %miu calculation
                    miu(l) = -1/C*v*rDiffMat(l,:)'/rDistances(l);
                    %frequency calculation
                    hzDopplerShift = Fc*miu(l);
                    %Al Matrix
                    %hermitianA(:,:,l) = diag(exp(2*pi*i*hzDopplerShift*(1:N)/Fs))';
                    %dopplerFixedSignal(l,:) = (hermitianA(:,:,l)*receivedSignal(l,:).').';
                    vHermitianA = exp(-2*pi*i*hzDopplerShift*(1:size(receivedSignal,2))/Fs);
                    dopplerFixedSignal(l,:) =(vHermitianA.*receivedSignal(l,:)).';
                end

                % Ntilda = N-sampleMaxDelay;
                %Evaluate V
                V = zeros(N,L);
                for l=1:L
                    V(:,l) = dopplerFixedSignal(l,(1+sampleDelay(l)):(N+sampleDelay(l)))';
                end

                %Evaluate The Cost
                if (method==2) %Unknown Signals
                    Qtilda = V'*V;
                    %Evaluate LambdaMax(Qtilda)
                    cost = max(eig(Qtilda));             
                elseif (method==1) %Known Signals
                    cost = originalSignal(1:N)*V*V'*originalSignal(1:N)';
                end
                if (cost>maxCost)
                    x = gX;
                    y = gY;
                    vx = gVx;
                    vy = gVy;
                    maxCost = cost;
                end
            end
        end
    end
end

close(waitbarHandle);
