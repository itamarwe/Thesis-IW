clear all;
close all;

%***** Init *****
Ns = 2^9;
Fs = 2^23; %Sampling Frequency [Hz]
Ts = 1/Fs;
T0=Ns/Fs; % observation time
time_vec = (0:Ns-1)*Ts; %Time Vector
freq_vec = (0:Ns-1)*(Fs/Ns);
f_delta = freq_vec(2)-freq_vec(1);

%**** Creating The Original Signal ****
NOF = 201; %Number of Frequencies
if NOF>1
    f0 = (1:NOF)./T0;
    f1 = f0(f0<0.4*Fs);
else
    f1 = 1/T0;
end

amplitudes = sin(30*pi*f1/f1(end))./(30*pi*f1/f1(end));

%amplitudes = randn(size(f1));

original_sig = amplitudes*cos( 2*pi*(f1')*(time_vec-T0/2)) ;

%******************************************************************
%**** Comparing the original signal with a shifted back signal ****
%******************************************************************

%Defining The Time and Frequency Shifts
DT = 1e-12; %[Sec]
DD = 100; %[Hz]

%Creating a signal shifted in DT and DD
shifted_sig = amplitudes*cos( 2*pi*(f1')*(time_vec-DT-T0/2)) ;
shifted_sig=(shifted_sig.').*exp(1i*2*pi*DD*time_vec');

%Shifting Back The shifted signal
shifted_back_sig = time_freq_shift(shifted_sig,Fs,-DT,-DD);

figure;
plot(time_vec*1e6,real(original_sig),'-b',time_vec*1e6, real(shifted_back_sig),'-r')
legend('Original Signal','Shifted by DT and DD and Shifted Back by DT and DD');
xlabel 'Time[\mu sec]';
grid
drawnow;
title 'The Original Signal Vs. A Shifted Back Signal';

%Calculating the noise induced by shifting one signal back and forth by
%various DT

DT = logspace(log10(T0*1e-6),log10(T0/2),20);
DD = 3000; %[Hz]
SNRdB = zeros(length(DT),1);
for DT_ind=1:length(DT)
    DT_ind
    %shifted_back_sig = time_freq_shift(time_freq_shift(original_sig,Fs,DT(DT_ind),DD),Fs,-DT(DT_ind),-DD);
    shifted_sig = amplitudes*cos( 2*pi*(f1')*(time_vec-DT(DT_ind)-T0/2)) ;
    shifted_sig=(shifted_sig.').*exp(1i*2*pi*DD*time_vec');
    shifted_back_sig = time_freq_shift(shifted_sig,Fs,-DT(DT_ind),-DD);
    SNRdB(DT_ind) = 10*log(norm(original_sig)/norm(shifted_back_sig-original_sig));
end%DD

figure;
semilogx(DT*1e6,SNRdB);
xlabel 'DT Time Shift [\mu sec]'
ylabel 'SNR[dB]'
grid
title 'Induced Noise Original Signal Vs. Shifted Back Signal';
drawnow;

%********************************************
%**** Comparing two shifted back signals ****
%********************************************

%Here the two signals are time and frequency shifted:
% - The first is shifted by DT0 and DD and then shifted back by DT0 and DD
% - The second is shifted by DT and DD and then shifted back by DT and DD

DT0 = T0/10;
DT = linspace(1e-6,T0/10,20);
DD = 3000; %[Hz]
SNRdB = zeros(length(DT),1);

%A signal shifted by DT0
shifted_sig1 = amplitudes*cos( 2*pi*(f1')*(time_vec-DT0-T0/2)) ;
shifted_sig1=(shifted_sig1.').*exp(1i*2*pi*DD*time_vec');
shifted_back_sig1= time_freq_shift(shifted_sig1,Fs,-DT0,-DD);
for DT_ind=1:length(DT)
    DT_ind
    %shifted_back_sig = time_freq_shift(time_freq_shift(original_sig,Fs,DT(DT_ind),DD),Fs,-DT(DT_ind),-DD);
    shifted_sig2 = amplitudes*cos( 2*pi*(f1')*(time_vec-DT(DT_ind)-DT0-T0/2)) ;
    shifted_sig2=(shifted_sig2.').*exp(1i*2*pi*DD*time_vec');
    shifted_back_sig2 = time_freq_shift(shifted_sig2,Fs,-DT(DT_ind)-DT0,-DD);
    SNRdB(DT_ind) = 10*log(norm(original_sig)/norm(shifted_back_sig2.'-shifted_back_sig1.'));
end%DD

figure;
plot(time_vec*1e6,real(shifted_back_sig1),'-b',time_vec*1e6, real(shifted_back_sig2),'-r')
legend('Shifted By DT0 and DD, and shifted back by DT0 and DD','Shifted By DT0+DT and DD, and shifted back by DT0+DT and DD');
xlabel 'Time[\mu sec]';
grid
drawnow;
title 'Shifting Two Signals Back And Forth';

figure;
semilogx(DT*1e6,SNRdB);
xlabel 'DT Time Shift[\mu sec]'
ylabel 'SNR[dB]'
grid
title 'Induced Noise Due Shifting Two Signals Back and Forth';
drawnow;

%**************************************************************
%**** Comparing a shifted signal and a shifted back signal ****
%**************************************************************

%Here the two signals are time and frequency shifted:
% - The first is shifted by DT0 and DD
% - The second is shifted by DT0+DT and DD and then shifted back by DT and DD

DT0 = T0/10;
DT = linspace(1e-6,T0/10,20);
DD = 3000; %[Hz]
SNRdB = zeros(length(DT),1);

%A signal shifted by DT0
shifted_sig1 = amplitudes*cos( 2*pi*(f1')*(time_vec-DT0-T0/2)) ;
shifted_sig1=(shifted_sig1.').*exp(1i*2*pi*DD*time_vec');
for DT_ind=1:length(DT)
    DT_ind
    %shifted_back_sig = time_freq_shift(time_freq_shift(original_sig,Fs,DT(DT_ind),DD),Fs,-DT(DT_ind),-DD);
    shifted_sig2 = amplitudes*cos( 2*pi*(f1')*(time_vec-DT(DT_ind)-DT0-T0/2)) ;
    shifted_sig2=(shifted_sig2.').*exp(1i*2*pi*DD*time_vec');
    shifted_back_sig2 = time_freq_shift(shifted_sig2,Fs,-DT(DT_ind),0);
    SNRdB(DT_ind) = 10*log(norm(original_sig)/norm(shifted_back_sig2.'-shifted_sig1));
end%DD

figure;
plot(time_vec*1e6,real(shifted_sig1),'-b',time_vec*1e6, real(shifted_back_sig2),'-r')
legend('Shifted by DT0 and DD','Shifted By DT0+DT and DD, and shifted back by DT');
xlabel 'Time[\mu sec]';
grid
drawnow;
title 'A Shifted Signal Vs. A Shifted Back Signal';

figure;
semilogx(DT*1e6,SNRdB);
xlabel 'DT Time Shift [\mu sec]'
ylabel 'SNR[dB]'
grid
title 'Induced Noise, Shifted Signal Vs. Shifted Back Signal';
drawnow;
