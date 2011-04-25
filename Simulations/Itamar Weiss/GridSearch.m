function cost = GridSearch(gridX, gridY, gridVx, gridVy,originalSignal, rReceiverMat, receivedSignal, Fc, C, Fs, N, method)
% cost = GridSearch(gridX, gridY, gridVx, gridVy, rTransmitter, v,rReceiverMat, Fc, C, Fs, method)
% method: 1 - UnknownSignals, 2 - KnownSignals, 3 - KnownNarrowBandSignals
%**********************************************************8
% GRID SEARCH *********************************8
% *************************************************8
TRIALS = length(gridX);
L = size(rReceiverMat,1);
lambdaMax = zeros(1,TRIALS);

Ntilde = size(receivedSignal,2);
k=0:Ntilde-1;
[kMat,nMat] = meshgrid(0:Ntilde-1,0:Ntilde-1);
%F = exp(-2*pi*i/Ntilde*kMat.*nMat);
%invF = 1/Ntilde*conj(F);
cost = zeros(1,TRIALS);

waitbarHandle = waitbar(0,'Processing');
for trial=1:TRIALS

    if mod(trial,10)==0
        waitbar(trial/TRIALS,waitbarHandle);
    end

    rTransmitter = [gridX(trial) gridY(trial)];
    v = [gridVx(trial) gridVy(trial)];

    rTransmitterMat = ones(L,2)*[rTransmitter(1) 0;0 rTransmitter(2)];
    rDiffMat = rTransmitterMat-rReceiverMat;
    rDistances = sqrt(rDiffMat(:,1).^2+rDiffMat(:,2).^2);
    %microsecDelay = rDistances/(C/1e6);
   % sampleDelay = microsecDelay*Fs/1e6;
    [maxDelay delay] = CalcMicroSecRxSignalDelay(rTransmitter, rReceiverMat,C);
    %Evaluate hermitian A
    %    hermitianA = zeros(N,N,L);
    %miu = zeros(L,1);
    dopplerFixedSignal = zeros(L,size(receivedSignal,2));
    for l=1:L
        %miu calculation
        miu = -1/C*v*rDiffMat(l,:)'/rDistances(l);
        %frequency calculation
        hzDopplerShift = Fc*miu;
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
        m = delay(l)*(Fs*1e-6); %04/07/10
        vk = fft(dopplerFixedSignal(l,:));     
       % vk = vk.*exp(2*pi*i/Ntilde*k*sampleDelay(l));% 04/07/10
        vk = vk.*exp(2*pi*i/Ntilde*k*m); % 04/07/10
        vl = ifft(vk);
        V(:,l) = vl(1:N);
        
        %vk = fft(dopplerFixedSignal(l,:));
        %vk = vk.*exp(2*pi*i/Ntilde*k*sampleDelay(l));
        %vl = ifft(vk);
        %V(:,l) = vl(1:N);
        
        %        V(:,l) = dopplerFixedSignal(l,(1+sampleDelay(l)):(N+sampleDelay(l)))';
        % V(:,l) = dopplerFixedSignal(l,(1+sampleMaxDelay-sampleDelay(l)):(N+sampleMaxDelay-sampleDelay(l)))';
    end

    %Evaluate Qtilda
    if (method==1)
        Qtilda = V'*V;
        %Evaluate LambdaMax(Qtilda)
        cost(trial) = max(eig(Qtilda));
    elseif (method==2)
        cost(trial) = originalSignal(1:N)*V*V'*originalSignal(1:N)';
    end
end
close(waitbarHandle);