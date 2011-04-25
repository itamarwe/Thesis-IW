function [DT DF DT0 DF0] = DT_DF_gridsearch(receivedSignal, rReceiverMat, rTransmitter,v,C,Fc,Fs)
%This Function get as an input the received signal in the receivers, and
%produces as an output the Differential Time, and Diffrential frequency
%between the receivers. Meaning: DT[i,j] = the estimated TOA diffrenece between the
%i-th receiver and the j-th receiver.

%The number of receivers
L = size(receivedSignal,1);
%The length of the signal [samples]
N = size(receivedSignal,2);
%The maximum possible value of the frequency shift
%maxDf = 2*Fc*norm(v)/C;
maxDf = 2*Fc*10000/C;

rTransmitterMat = ones(L,2)*[rTransmitter(1) 0;0 rTransmitter(2)];
rDiffMat = rTransmitterMat-rReceiverMat;

%The differential distance matrix
rDistances = sqrt(rDiffMat(:,1).^2+rDiffMat(:,2).^2);
for i=1:L
    for j=(i+1):L
        rDifDistances(i,j) = abs(rDistances(i)-rDistances(j));
    end
end

%The Differential delay matrix [samples]
sampleDelay = (rDifDistances/(C/Fs));
%The maximum differential delay [samples]
%maxDt = max(max(sampleDelay));
maxDt = 2000/(C/Fs);
%DT0 - The DT matrix if there was no noise, calculated accroding to the
%geometry.
%DF0 - The DF matrix if there was no noise, calculated accroding to the
%geometry.
DT0 =zeros(L);
DF0 =zeros(L);
DT =zeros(L);
DF =zeros(L);

k_vec = 0:N-1;

for k=1:1
    for l = (k+1):L
        %DT0 calculation
        DT0(k,l) = norm((rReceiverMat(k,:)-rTransmitter))/C-norm((rReceiverMat(l,:)-rTransmitter))/C;
        %Because DT is anti-symmetric
        %%DT0(l,k) =    - DT0(k,l);
        %DF0 calculation
        DF0(k,l) = -Fc/C*((v*(rTransmitter-rReceiverMat(k,:))')/norm((rTransmitter-rReceiverMat(k,:))))+Fc/C*((v*(rTransmitter-rReceiverMat(l,:))')/norm((rTransmitter-rReceiverMat(l,:))));
        %Because DF is anti-symmetric
        %%DF0(l,k) =    - DF0(k,l);
    end
end


for k=1:1
    for l=(k+1):L
        dt = DT0(k,l)*Fs;
        df = DF0(k,l);
        dtScale = maxDt*2;
        dfScale = maxDf;

        dt_vec = linspace(-dtScale/2,dtScale/2,100)+dt+(rand-0.5)*1e-4;
        df_vec = linspace(-dfScale/2,dfScale/2,100)+df+(rand-0.5)*5;
    
        flag=0;
        while dtScale>(1e-4/C*Fs)
            A_mat = exp(-2*pi*1i*(df_vec'*(0:(N-1)))/Fs);
            CAF = zeros(length(dt_vec),length(df_vec));
            for i=1:length(dt_vec)
                %Time shift of the k-th signal
                %ak = fft(receivedSignal(k,:));
                %bk = ak.*exp(2*pi*1i*(k_vec/N)*dt_vec(i));
                %delayedSignal = ifft(bk);
                delayedSignal = time_freq_shift(receivedSignal(k,:).',Fs,-dt_vec(i)/Fs,0);
                %df = df_vec(2)-df_vec(1);
                %A_vec0 = exp(2*pi*1i*df_vec(1)*(0:(N-1))/Fs);
                %A_vec = exp(2*pi*1i*df*(0:(N-1))/Fs);
                %Give an intial freq shift to the delayed signal
                % shiftedSignal = A_vec0.*delayedSignal;
                for j= 1:length(df_vec)
                    %Freq shift of the k-th signal
                    %A_vec = exp(-2*pi*1i*df_vec(j)*((0:(N-1))/Fs));
                    %shiftedSignal = A_vec.*delayedSignal;
                    %%shiftedSignal = A_mat(j,:).*delayedSignal;
                    shiftedSignal = time_freq_shift(delayedSignal,Fs,0,-df_vec(j));
                    CAF(i,j) = shiftedSignal*receivedSignal(l,:)';
                    %shiftedSignal = A_vec.*shiftedSignal;
                    %    plot(abs(shiftedSignal),'r');hold on; plot(abs(receivedSignal(l,:)));hold off;
                end
            end
            
            if(flag)
            surf(dt_vec,df_vec,abs(CAF));
             title(['k=' num2str(k) ', l=' num2str(l)])
             keyboard
             flag=0;
            end
            [indt,indf] = find(abs(CAF) == max(max(abs(CAF))));

            df = df_vec(indf(1));
            dt = dt_vec(indt(1));
            dtScale = dtScale*0.6;
            dfScale = dfScale*0.6;
            dt_vec = [linspace(-dtScale/2,dtScale/2,3)]+dt;
            df_vec = [linspace(-dfScale/2,dfScale/2,3)]+df;

        end

           
        if (indf(1)~=1)&&(indf(1)~=length(df_vec))
            df = peak_identification(df_vec,abs(CAF(indt(1),:)),indf(1));
        end
        if (indt(1)~=1)&&(indt(1)~=length(dt_vec))
            dt = peak_identification(dt_vec,abs(CAF(:,indf(1))),indt(1));
        end

        DF(k,l) = df;

        DT(k,l) = dt/Fs;
    end
end

