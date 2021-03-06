function [x y vx vy] = GridSearch5(gridX, gridY, gridVx, gridVy,originalSignal, rReceiverMat, receivedSignal, Fc, C, Fs, N, method,DT,DF,CRB_toa,CRB_foa)
% cost = GridSearch(gridX, gridY, gridVx, gridVy, rTransmitter, v,rReceiverMat, Fc, C, Fs, method)
% method: 1 - UnknownSignals, 2 - KnownSignals, 3 - DPD no doppler Unknown
% 4 - DPD no doppler Known 5 - Coventional FDOA-TDOA unknown
%**********************************************************8
% GRID SEARCH *********************************8
% *************************************************8

%The number of receivers
L = size(rReceiverMat,1);

%The maximum cost calculated in the currrent grid search
maxCost = -inf;

%The number of samples in the received signal
Ntilde = size(receivedSignal,2);

%A vector used for time translations in the frequency domain
k=0:Ntilde-1;

if (method==5)
    %The weights for the cost function
    w_toa = 1/CRB_toa;
    w_foa = 1/CRB_foa;
end

for gX = gridX
    for gY = gridY
        for gVx = gridVx
            for gVy = gridVy

                %The examined position of the transmitter
                rTransmitter = [gX gY];
                %The examined velocity of the transmitter
                v = [gVx gVy];
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

                        if (method==3)||(method==4)
                            %For DPD No doppler Method
                            dopplerFixedSignal(l,:) =receivedSignal(l,:).';
                        end
                    end
                    %Evaluate V
                    V = zeros(N,L);
                    for l=1:L
                        m = delay(l)*(Fs*1e-6); %04/07/10
                        vk = fft(dopplerFixedSignal(l,:));
                        % vk = vk.*exp(2*pi*i/Ntilde*k*sampleDelay(l));% 04/07/10
                        vk = vk.*exp(2*pi*i/Ntilde*k*m); % 04/07/10
                        vl = ifft(vk);
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
                        cost = originalSignal(1:N)*V*V'*originalSignal(1:N)';
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
                        %The cost functino of the conventional method
                        %The minus sign is put here so that maximum of the cost
                        %function (with the minus) will be in the minimum of
                        %the cost function (without the minus)
                        cost = -(w_toa*sum(sum((DT-DT0).^2))+w_foa*sum(sum((DF-DF0).^2)));
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