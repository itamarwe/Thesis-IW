clear all;
close all;


% ******Initialize the scenario
disp('Initializing the scenario');

C = 3e08; %Speed of Signal Propagation (Usually speed of light)[m/s]
Ns = 512; %Number of samples
Fs = 2^23; %Sampling Frequency [Hz]
Ts = 1/Fs;
T0=Ns/Fs; % observation time
L = 6; %Number of receivers
R = 1000;% [m] Receivers Distance From the axis center
W = 30; % [samples] Signal width

time_vec = (0:Ns-1)*Ts; %Time Vector
freq_vec = (0:Ns-1)*(Fs/Ns);
f_delta = freq_vec(2)-freq_vec(1);

rTransmitter = [100,300]; %Transmitter position [x,y] [m]
v = [200,200];%Transmitter velocity [vx, vy] [m/s]
Fc = 1e9; % Carrier Frequency[Hz]
B = 0.4*Fs; %Transmitter Bandwidth
receiverGeometry = 1; %1 - Circular 2- Linear

SNRdB = 100;


%Receivers Positions
rReceiverMat = LocateReceivers(R,L,receiverGeometry);

%Plot the scenario geometry
figure
plot(rReceiverMat(:,1),rReceiverMat(:,2),'k*',...
    rTransmitter(1),rTransmitter(2),'kp') ; grid on
xlabel('[m]');
ylabel('[m]')
axis equal
legend('Receivers','Transmitter');
title 'Scenario Geometry';
drawnow

% ******Create Shifted Signal
disp('Creating the original signal');
%Create Original Signal
NOF = 201; %Number of Frequencies
if NOF>1
    f0 = (1:NOF)./T0;
    f1 = f0(f0<B);
else
    f1 = 1/T0;
end

amplitudes = sin(W*pi*f1/f1(end))./(W*pi*f1/f1(end));
%amplitudes = randn(size(f1));
int = amplitudes*cos( 2*pi*(f1')*(time_vec-T0/2)) ;

figure
plot(time_vec*1e6,int)
xlabel('Time [\mu Sec]')
ylabel('Amplitude')
title('Signal Waveform')
grid
drawnow

%Calculate shifts
rTransmitterMat = ones(L,2)*[rTransmitter(1) 0;0 rTransmitter(2)];
rDiffMat = rTransmitterMat-rReceiverMat; %The difference matrix
rDistances = sqrt(rDiffMat(:,1).^2+rDiffMat(:,2).^2); %The Distances Vector

TOA =  rDistances/C;
Doppler_shift = zeros(L,1);
for ell=1:L
    Doppler_shift(ell) = -Fc/C*v*rDiffMat(ell,:)'/rDistances(ell);
end;

Delayed_sig = zeros(Ns,L);
Delayed_sig0 = Delayed_sig;
s_dot = zeros(Ns,L);
sig = amplitudes*cos( 2*pi*(f1')*(time_vec-T0/2)) ;
for ell = 1:L
    DT = TOA(ell) - TOA(1);
    DD = Doppler_shift(ell) - Doppler_shift(1);
    int(ell,:) = amplitudes*cos( 2*pi*(f1')*(time_vec-DT-T0/2)) ;
    int_dot(ell,:) = -amplitudes*2*pi*diag(f1)*sin( 2*pi*(f1')*(time_vec-DT-T0/2));
    int(ell,:)=int(ell,:)./norm(int(ell,:));
    int_dot(ell,:) = int_dot(ell,:)./norm(int(ell,:));
    int0(ell,:) = amplitudes*cos( 2*pi*(f1')*(time_vec-TOA(ell)-T0/2)) ;
    int0_dot(ell,:) = -amplitudes*2*pi*diag(f1)*sin( 2*pi*(f1')*(time_vec-TOA(ell)-T0/2));
    int0(ell,:) = int0(ell,:)./norm(int0(ell,:));
    int0_dot(ell,:) = int0_dot(ell,:)./norm(int0(ell,:));
    Delayed_sig(:,ell)=(int(ell,:).').*exp(1i*2*pi*DD*time_vec');
    Delayed_sig0(:,ell)=(int0(ell,:).').*exp(1i*2*pi*Doppler_shift(ell)*time_vec');

end;%ell


%Add Noise
sigma2Vec = 10.^(-SNRdB/10);% noise variance

disp('Calculating the cost functions');
r_vec_t = zeros(Ns,L);
r_vec_t0 = zeros(Ns,L);
for snr_ind = 1:length(SNRdB)
    tic
    %Add Noise to the signals
    sigma2  = sigma2Vec(snr_ind);
    sigma_sq = sigma2/Ns;
    for ell = 1:L
        noise = sqrt(sigma_sq/2)*(randn(Ns,1)+1i*randn(Ns,1));
        r_vec_t(:,ell) =  Delayed_sig(:,ell) + noise ;
    end;%ell

    %Calculate CRB for the conventional method
    [CRB_toa, CRB_foa] = CRB_toa_foa(Fc,C, int0, int0_dot, Ns, Fs, ones(1,L), rTransmitter,v, rReceiverMat,SNRdB(snr_ind));
    CRB_toa = max(CRB_toa);
    CRB_foa = max(CRB_foa);
    w_toa = 1/CRB_toa; %TOA Measurements weight for conventional method
    w_foa = 1/CRB_foa; %FOA Measurements weight for conventional method

    [DT_conv,DF_conv] = DT_DF_gridsearch(r_vec_t.', rReceiverMat, rTransmitter,v,C,Fc,Fs);

    %****** Calculate Cost Function
    scale = .01;%[m]
    Xsearch = linspace(-scale/2, scale/2,50)+rTransmitter(1);
    Ysearch = linspace(-scale/2,scale/2,50)+rTransmitter(2);
    CF = zeros(length(Xsearch),length(Ysearch)); %The Unknown Signals Cost Function
    CF_known = zeros(length(Xsearch),length(Ysearch)); %The Known Signals Cost Function
    CF_conv = zeros(length(Xsearch),length(Ysearch)); %The Conventional Unknown Signals Cost Function
    for x_ind = 1:length(Xsearch)
        for y_ind = 1:length(Ysearch)
            %Calculate shifts
            rTransmitterMat = ones(L,2)*[Xsearch(x_ind) 0;0 Ysearch(y_ind)];
            rDiffMat = rTransmitterMat-rReceiverMat; %The difference matrix
            rDistances = sqrt(rDiffMat(:,1).^2+rDiffMat(:,2).^2); %The Distances Vector

            TOA =  rDistances/C;
            Doppler_shift = zeros(L,1);
            for ell=1:L
                Doppler_shift(ell) = -Fc/C*v*rDiffMat(ell,:)'/rDistances(ell);
            end;

            Vk=zeros(Ns,L);
            for ell = 1:L

                DD = Doppler_shift(ell) - Doppler_shift(1);
                DT = TOA(ell)-TOA(1);
                Vk(:,ell) = time_freq_shift(r_vec_t(:,ell),Fs,-DT,-DD);

            end;%ell
            Lambda  = svd(Vk'*Vk); % (Eq.22 - old)
            CF_temp = abs(Lambda(1));
            CF(x_ind,y_ind) = CF_temp; %Unknown Signals Cost Function
            CF_known(x_ind,y_ind) = sig*Vk*(Vk'*sig');%Known Signals Cost function

            %Conventional method
            DT0 = zeros(L);
            DF0 = zeros(L);
            for k2=1:L
                for l2 = (k2+1):L
                    %DT0 calculation
                    %DT0(k2,l2) = norm((rReceiverMat(k2,:)-rTransmitter))/C-norm((rReceiverMat(l2,:)-rTransmitter))/C;
                    DT0(k2,l2) = (rDistances(k2)-rDistances(l2))/C;
                    %DF0 calculation
                    %DF0(k2,l2) = -Fc/C*((v*(rTransmitter-rReceiverMat(k2,:))')/norm((rTransmitter-rReceiverMat(k2,:))))+Fc/C*((v*(rTransmitter-rReceiverMat(l2,:))')/norm((rTransmitter-rReceiverMat(l2,:))));
                    DF0(k2,l2) = Doppler_shift(k2)-Doppler_shift(l2);
                end
            end
            %The cost function of the conventional method
            %The minus sign is put here so that maximum of the cost
            %function (with the minus) will be in the minimum of
            %the cost function (without the minus)
            CF_conv(x_ind,y_ind) = -(w_toa*sum(sum((DT_conv(1,:)-DT0(1,:)).^2))+w_foa*sum(sum((DF_conv(1,:)-DF0(1,:)).^2)));
        end; % y_ind
    end; % x_ind

    figure
    surf(Xsearch, Ysearch,CF), shading interp
    title('Unknown DPD CF')
    xlabel('[m]');
    ylabel('[m]')

    figure
    contour(Xsearch, Ysearch, CF),
    title('Unknown DPD CF')
    xlabel('[m]');
    ylabel('[m]')

    figure
    surf(Xsearch, Ysearch,CF_known), shading interp
    title('Known DPD CF')
    xlabel('[m]');
    ylabel('[m]')

    figure
    contour(Xsearch, Ysearch, CF_known),
    title('Known DPD CF')
    xlabel('[m]');
    ylabel('[m]')

    figure
    surf(Xsearch, Ysearch,CF_conv), shading interp
    title('Conventional CF')
    xlabel('[m]');
    ylabel('[m]')

    figure
    contour(Xsearch, Ysearch, CF_conv),
    title('Conventional CF')
    xlabel('[m]');
    ylabel('[m]')

    drawnow;

    %Calculate Velocity Cost Function
    scale = .01;%[m/s]
    Vxsearch = linspace(-scale/2, scale/2,50)+v(1);
    Vysearch = linspace(-scale/2,scale/2,50)+v(2);
    CF = zeros(length(Vxsearch),length(Vysearch)); %The Cost Function
    for Vx_ind = 1:length(Vxsearch)
        for Vy_ind = 1:length(Vysearch)
            %Calculate shifts
            rTransmitterMat = ones(L,2)*[rTransmitter(1) 0;0 rTransmitter(2)];
            rDiffMat = rTransmitterMat-rReceiverMat; %The difference matrix
            rDistances = sqrt(rDiffMat(:,1).^2+rDiffMat(:,2).^2); %The Distances Vector

            TOA =  rDistances/C;
            Doppler_shift = zeros(L,1);
            for ell=1:L
                Doppler_shift(ell) = -Fc/C*[Vxsearch(Vx_ind),Vysearch(Vy_ind)]*rDiffMat(ell,:)'/rDistances(ell);
            end;

            Vk=zeros(Ns,L);
            for ell = 1:L

                DD = Doppler_shift(ell) - Doppler_shift(1);
                DT = TOA(ell)-TOA(1);
                Vk(:,ell) = time_freq_shift(r_vec_t(:,ell),Fs,-DT,-DD);

            end;%ell
            Lambda  = svd(Vk'*Vk); % (Eq.22 - old)
            CF_temp = abs(Lambda(1));
            CF(Vx_ind,Vy_ind) = CF_temp;
            CF_known(Vx_ind,Vy_ind) = sig*Vk*(Vk'*sig');%Known Signals Cost function
            %Conventional method
            DT0 = zeros(L);
            DF0 = zeros(L);
            for k2=1:L
                for l2 = (k2+1):L
                    %DT0 calculation
                    %DT0(k2,l2) = norm((rReceiverMat(k2,:)-rTransmitter))/C-norm((rReceiverMat(l2,:)-rTransmitter))/C;
                    DT0(k2,l2) = (rDistances(k2)-rDistances(l2))/C;
                    %DF0 calculation
                    %DF0(k2,l2) = -Fc/C*((v*(rTransmitter-rReceiverMat(k2,:))')/norm((rTransmitter-rReceiverMat(k2,:))))+Fc/C*((v*(rTransmitter-rReceiverMat(l2,:))')/norm((rTransmitter-rReceiverMat(l2,:))));
                    DF0(k2,l2) = Doppler_shift(k2)-Doppler_shift(l2);
                end
            end
            %The cost function of the conventional method
            %The minus sign is put here so that maximum of the cost
            %function (with the minus) will be in the minimum of
            %the cost function (without the minus)
            CF_conv(Vx_ind,Vy_ind) = -(w_toa*sum(sum((DT_conv(1,:)-DT0(1,:)).^2))+w_foa*sum(sum((DF_conv(1,:)-DF0(1,:)).^2)));
        end; % y_ind
    end; % x_ind

    figure
    surf(Vxsearch, Vysearch,CF), shading interp
    title('Unknown  DPD Velocity CF')
    xlabel('[m/s]');
    ylabel('[m/s]')

    figure
    contour(Vxsearch, Vysearch, CF),
    title('Unknown DPD Velocity CF')
    xlabel('[m/s]');
    ylabel('[m/s]')
    
    figure
    surf(Vxsearch, Vysearch,CF_known), shading interp
    title('Known DPD Velocity CF')
    xlabel('[m/s]');
    ylabel('[m/s]')
    figure
    
    contour(Vxsearch, Vysearch, CF_known),
    title('Known DPD Velocity CF')
    xlabel('[m/s]');
    ylabel('[m/s]')
    
    figure
    surf(Vxsearch, Vysearch,CF_conv), shading interp
    title('Conventional Velocity CF')
    xlabel('[m/s]');
    ylabel('[m/s]')
    
    figure
    contour(Vxsearch, Vysearch, CF_conv),
    title('Conventional Velocity CF')
    xlabel('[m/s]');
    ylabel('[m/s]')
    drawnow

end

keyboard
close all;
drawnow;

%******************************
%**** Perform Grid Search *****
%******************************
SNRdB = 25:-10:-25;
ESTIMATIONS = 10;
PSCALE = 5000;
PRESOLUTION = 0.1;
VSCALE = 10000;


sigma2Vec = 10.^(-SNRdB/10);% noise variance

disp('Starting Grid Search');

for snr_ind = 1:length(SNRdB)
    sigma2  = sigma2Vec(snr_ind);
    sigma_sq = sigma2/Ns;
    estimationResultsKnown = zeros(ESTIMATIONS,4);
    estimationResultsUnknown = zeros(ESTIMATIONS,4);
    estimationResultsConv = zeros(ESTIMATIONS,4);

    for estimation = 1:ESTIMATIONS
        tic
        disp(['Experiment No. ' num2str(estimation) '/' num2str(ESTIMATIONS)]);
        disp(['SNRdB: ' num2str(SNRdB(snr_ind))]);

        %Add Noise to the signals
        for ell = 1:L
            noise = sqrt(sigma_sq/2)*(randn(Ns,1)+1i*randn(Ns,1));
            r_vec_t(:,ell) =  Delayed_sig(:,ell) + noise ;
        end;%ell

        %Calculate CRB for the conventional method
        [CRB_toa, CRB_foa] = CRB_toa_foa(Fc,C, int0, int0_dot, Ns, Fs, ones(1,L), rTransmitter,v, rReceiverMat,SNRdB(snr_ind));
        CRB_toa = max(CRB_toa);
        CRB_foa = max(CRB_foa);
        w_toa = 1/CRB_toa; %TOA Measurements weight for conventional method
        w_foa = 1/CRB_foa; %FOA Measurements weight for conventional method

        disp('Performing DT, DF Search...')
        [DT_conv,DF_conv] = DT_DF_gridsearch(r_vec_t.', rReceiverMat, rTransmitter,v,C,Fc,Fs);

        disp('Performing Unknown DPD Search...')
        pScale =  PSCALE;
        vScale = VSCALE;
        gridX = linspace(-pScale/2, pScale/2,10)+rTransmitter(1)+ (rand-.5)*pScale/5;
        gridY = linspace(-pScale/2, pScale/2,10)+rTransmitter(2)+(rand-.5)*pScale/5;
        gridVx = linspace(-vScale/2, vScale/2,10)+v(1)+(rand-.5)*vScale/5;
        gridVy = linspace(-vScale/2, vScale/2,10)+v(2)+(rand-.5)*vScale/5;

        while pScale > PRESOLUTION

            maxCost = -inf;

            for gX = gridX
                for gY = gridY
                    for gVx = gridVx
                        for gVy = gridVy

                            rTransmitterMat = ones(L,2)*[gX 0;0 gY];
                            vg = [gVx, gVy];
                            rDiffMat = rTransmitterMat-rReceiverMat; %The difference matrix
                            rDistances = sqrt(rDiffMat(:,1).^2+rDiffMat(:,2).^2); %The Distances Vector

                            TOA =  rDistances/C;
                            Doppler_shift = zeros(L,1);
                            for ell=1:L
                                Doppler_shift(ell) = -Fc/C*vg*rDiffMat(ell,:)'/rDistances(ell);
                            end;

                            Vk=zeros(Ns,L);
                            for ell = 1:L

                                DD = Doppler_shift(ell) - Doppler_shift(1);
                                DT = TOA(ell)-TOA(1);
                                Vk(:,ell) = time_freq_shift(r_vec_t(:,ell),Fs,-DT,-DD);

                            end;%ell
                            Lambda  = svd(Vk'*Vk); % (Eq.22 - old)
                            cost = abs(Lambda(1));

                            if (cost>maxCost)
                                x=gX;
                                y=gY;
                                vx=gVx;
                                vy=gVy;
                                maxCost = cost;
                            end
                        end%gVy
                    end%gVx
                end%gY
            end%gX

            pScale=pScale*.6;
            vScale=vScale*.6;
            gridX = [linspace(-pScale/2, pScale/2,3)] + x;
            gridY = [linspace(-pScale/2, pScale/2,3)] + y;
            gridVx = [linspace(-vScale/2, vScale/2,3)]+vx;
            gridVy = [linspace(-vScale/2, vScale/2,3)]+vy;
        end
        disp('Unknown DPD estimation:');
        disp(['X:' num2str(x) ' Y:' num2str(y) ' Vx:' num2str(vx) ' Vy:' num2str(vy)]);

        estimationResultsUnknown(estimation,:) = [x y vx vy];
        estimationResultsUnknown(estimation,:) = estimationResultsUnknown(estimation,:)-[rTransmitter v];

        disp('Performing Known DPD Search...')
        pScale =  PSCALE;
        vScale = VSCALE;
        gridX = linspace(-pScale/2, pScale/2,10)+rTransmitter(1)+ (rand-.5)*pScale/5;
        gridY = linspace(-pScale/2, pScale/2,10)+rTransmitter(2)+(rand-.5)*pScale/5;
        gridVx = linspace(-vScale/2, vScale/2,10)+v(1)+(rand-.5)*vScale/5;
        gridVy = linspace(-vScale/2, vScale/2,10)+v(2)+(rand-.5)*vScale/5;

        while pScale > PRESOLUTION

            maxCost = -inf;

            for gX = gridX
                for gY = gridY
                    for gVx = gridVx
                        for gVy = gridVy

                            rTransmitterMat = ones(L,2)*[gX 0;0 gY];
                            vg = [gVx, gVy];
                            rDiffMat = rTransmitterMat-rReceiverMat; %The difference matrix
                            rDistances = sqrt(rDiffMat(:,1).^2+rDiffMat(:,2).^2); %The Distances Vector

                            TOA =  rDistances/C;
                            Doppler_shift = zeros(L,1);
                            for ell=1:L
                                Doppler_shift(ell) = -Fc/C*vg*rDiffMat(ell,:)'/rDistances(ell);
                            end;

                            Vk=zeros(Ns,L);
                            for ell = 1:L

                                DD = Doppler_shift(ell) - Doppler_shift(1);
                                DT = TOA(ell)-TOA(1);
                                Vk(:,ell) = time_freq_shift(r_vec_t(:,ell),Fs,-DT,-DD);

                            end;%ell

                            cost = sig*Vk*(Vk'*sig');%Known Signals Cost function
                            if (cost>maxCost)
                                x=gX;
                                y=gY;
                                vx=gVx;
                                vy=gVy;
                                maxCost = cost;
                            end
                        end%gVy
                    end%gVx
                end%gY
            end%gX

            pScale=pScale*.6;
            vScale=vScale*.6;
            gridX = [linspace(-pScale/2, pScale/2,3)] + x;
            gridY = [linspace(-pScale/2, pScale/2,3)] + y;
            gridVx = [linspace(-vScale/2, vScale/2,3)]+vx;
            gridVy = [linspace(-vScale/2, vScale/2,3)]+vy;
        end
        disp('Known DPD estimation:');
        disp(['X:' num2str(x) ' Y:' num2str(y) ' Vx:' num2str(vx) ' Vy:' num2str(vy)]);

        estimationResultsKnown(estimation,:) = [x y vx vy];
        estimationResultsKnown(estimation,:) = estimationResultsKnown(estimation,:)-[rTransmitter v];

        disp('Performing Conventional Search...')
        pScale =  PSCALE;
        vScale = VSCALE;
        gridX = linspace(-pScale/2, pScale/2,10)+rTransmitter(1)+ (rand-.5)*pScale/5;
        gridY = linspace(-pScale/2, pScale/2,10)+rTransmitter(2)+(rand-.5)*pScale/5;
        gridVx = linspace(-vScale/2, vScale/2,10)+v(1)+(rand-.5)*vScale/5;
        gridVy = linspace(-vScale/2, vScale/2,10)+v(2)+(rand-.5)*vScale/5;

        while pScale > PRESOLUTION
            maxCost = -inf;

            for gX = gridX
                for gY = gridY
                    for gVx = gridVx
                        for gVy = gridVy

                            rTransmitterMat = ones(L,2)*[gX 0;0 gY];
                            vg = [gVx, gVy];
                            rDiffMat = rTransmitterMat-rReceiverMat; %The difference matrix
                            rDistances = sqrt(rDiffMat(:,1).^2+rDiffMat(:,2).^2); %The Distances Vector

                            TOA =  rDistances/C;
                            Doppler_shift = zeros(L,1);
                            for ell=1:L
                                Doppler_shift(ell) = -Fc/C*vg*rDiffMat(ell,:)'/rDistances(ell);
                            end;



                            DT0 = zeros(L);
                            DF0 = zeros(L);
                            for k2=1:L
                                for l2 = (k2+1):L
                                    %DT0 calculation
                                    %DT0(k2,l2) = norm((rReceiverMat(k2,:)-rTransmitter))/C-norm((rReceiverMat(l2,:)-rTransmitter))/C;
                                    DT0(k2,l2) = (rDistances(k2)-rDistances(l2))/C;
                                    %DF0 calculation
                                    %DF0(k2,l2) = -Fc/C*((v*(rTransmitter-rReceiverMat(k2,:))')/norm((rTransmitter-rReceiverMat(k2,:))))+Fc/C*((v*(rTransmitter-rReceiverMat(l2,:))')/norm((rTransmitter-rReceiverMat(l2,:))));
                                    DF0(k2,l2) = Doppler_shift(k2)-Doppler_shift(l2);
                                end
                            end
                            %The cost function of the conventional method
                            %The minus sign is put here so that maximum of the cost
                            %function (with the minus) will be in the minimum of
                            %the cost function (without the minus)
                            cost = -(w_toa*sum(sum((DT_conv(1,:)-DT0(1,:)).^2))+w_foa*sum(sum((DF_conv(1,:)-DF0(1,:)).^2)));

                            if (cost>maxCost)
                                x=gX;
                                y=gY;
                                vx=gVx;
                                vy=gVy;
                                maxCost = cost;
                            end
                        end%gVy
                    end%gVx
                end%gY
            end%gX

            pScale=pScale*.6;
            vScale=vScale*.6;
            gridX = [linspace(-pScale/2, pScale/2,3)] + x;
            gridY = [linspace(-pScale/2, pScale/2,3)] + y;
            gridVx = [linspace(-vScale/2, vScale/2,3)]+vx;
            gridVy = [linspace(-vScale/2, vScale/2,3)]+vy;
        end
        disp('Conventional estimation:');
        disp(['X:' num2str(x) ' Y:' num2str(y) ' Vx:' num2str(vx) ' Vy:' num2str(vy)]);

        estimationResultsConv(estimation,:) = [x y vx vy];
        estimationResultsConv(estimation,:) = estimationResultsConv(estimation,:)-[rTransmitter v];
        toc
    end%estimation
    scenarioResultsUnknown(snr_ind,1:4) = sqrt(mean(estimationResultsUnknown.^2));
    scenarioResultsUnknown(snr_ind,5) = sqrt(scenarioResultsUnknown(snr_ind,1)^2+scenarioResultsUnknown(snr_ind,2)^2);
    scenarioResultsUnknown(snr_ind,6) = sqrt(scenarioResultsUnknown(snr_ind,3)^2+scenarioResultsUnknown(snr_ind,4)^2);

    scenarioResultsKnown(snr_ind,1:4) = sqrt(mean(estimationResultsKnown.^2));
    scenarioResultsKnown(snr_ind,5) = sqrt(scenarioResultsKnown(snr_ind,1)^2+scenarioResultsKnown(snr_ind,2)^2);
    scenarioResultsKnown(snr_ind,6) = sqrt(scenarioResultsKnown(snr_ind,3)^2+scenarioResultsKnown(snr_ind,4)^2);

    scenarioResultsConv(snr_ind,1:4) = sqrt(mean(estimationResultsConv.^2));
    scenarioResultsConv(snr_ind,5) = sqrt(scenarioResultsConv(snr_ind,1)^2+scenarioResultsConv(snr_ind,2)^2);
    scenarioResultsConv(snr_ind,6) = sqrt(scenarioResultsConv(snr_ind,3)^2+scenarioResultsConv(snr_ind,4)^2);

    save experiment scenarioResultsKnown scenarioResultsUnknown scenarioResultsConv SNRdB;
    disp(' ');
end%SNRdB

semilogy(SNRdB,scenarioResultsUnknown(:,5),'b-s',SNRdB,scenarioResultsKnown(:,5),'r-o',SNRdB,scenarioResultsConv(:,5),'b-p');
xtext='SNR [dB]';
xlabel(xtext,'interpreter','latex','fontsize',16);
ytext='RMSE [m]';
ylabel(ytext,'interpreter','latex','fontsize',16);
grid;
set(gca,'xtick',-25:5:25);
set(gca,'xticklabel',-25:5:25);
legend('DPD Unknown Signals','DPD known Signals',...
    'Conventional Unknown Signals');
title('Position Estimation Preformance - Circular Array, Pulse Signal','interpreter','latex','fontsize',16);

figure;
semilogy(SNRdB,scenarioResultsUnknown(:,6),'b-s',SNRdB,scenarioResultsKnown(:,6),'r-o',SNRdB,scenarioResultsConv(:,6),'b-p');
xtext='SNR [dB]';
xlabel(xtext,'interpreter','latex','fontsize',16);
ytext='RMSE [m/s]';
ylabel(ytext,'interpreter','latex','fontsize',16);
grid;
set(gca,'xtick',-25:5:25);
set(gca,'xticklabel',-25:5:25);
legend('DPD Unknown Signals','DPD known Signals',...
    'Conventional Unknown Signals');
title('Velocity Estimation Preformance - Circular Array, Pulse Signal','interpreter','latex','fontsize',16);