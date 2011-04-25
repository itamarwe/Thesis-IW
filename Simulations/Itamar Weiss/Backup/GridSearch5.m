function [x y vx vy] = GridSearch5(gridX, gridY, gridVx, gridVy,originalSignal, rReceiverMat, receivedSignal, Fc, C, Fs, N, method)
% cost = GridSearch(gridX, gridY, gridVx, gridVy, rTransmitter, v,rReceiverMat, Fc, C, Fs, method)
% method: 1 - UnknownSignals, 2 - KnownSignals, 3 - KnownNarrowBandSignals
%**********************************************************8
% GRID SEARCH *********************************8
% *************************************************8
%TRIALS = length(gridX);
L = size(rReceiverMat,1);
%lambdaMax = zeros(1,TRIALS);
maxCost = -inf;

Ntilde = size(receivedSignal,2);
k=0:Ntilde-1;
%[kMat,nMat] = meshgrid(0:Ntilde-1,0:Ntilde-1);
%F = exp(-2*pi*i/Ntilde*kMat.*nMat);
%invF = 1/Ntilde*conj(F);

for gX = gridX
    for gY = gridY
        for gVx = gridVx
            for gVy = gridVy
                


    rTransmitter = [gX gY];
    v = [gVx gVy];

    rTransmitterMat = ones(L,2)*[rTransmitter(1) 0;0 rTransmitter(2)];
    rDiffMat = rTransmitterMat-rReceiverMat;
    rDistances = sqrt(rDiffMat(:,1).^2+rDiffMat(:,2).^2);
    microsecDelay = rDistances/(C/1e6);
    sampleDelay = microsecDelay*Fs/1e6;
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
        vk = fft(dopplerFixedSignal(l,:));
        vk = vk.*exp(2*pi*i/Ntilde*k*sampleDelay(l));
        vl = ifft(vk);
        V(:,l) = vl(1:N);
        %        V(:,l) = dopplerFixedSignal(l,(1+sampleDelay(l)):(N+sampleDelay(l)))';
        % V(:,l) = dopplerFixedSignal(l,(1+sampleMaxDelay-sampleDelay(l)):(N+sampleMaxDelay-sampleDelay(l)))';
    end

    %Evaluate Qtilda
    if (method==1)
        Qtilda = V'*V;
        %Evaluate LambdaMax(Qtilda)
        cost = max(eig(Qtilda));
    elseif (method==2)
        cost = originalSignal(1:N)*V*V'*originalSignal(1:N)';
    end
    
    if (cost>maxCost)
        x=gX;
        y=gY;
        vx=gVx;
        vy=gVy;
        maxCost = cost;
    end
            end %gVy
        end%gVx
    end%gY
end%gX