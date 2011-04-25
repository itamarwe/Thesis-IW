clear all;
close all;

Ns = 2^9;
Fs = 2^23; %Sampling Frequency [Hz]
Ts = 1/Fs;
T0=Ns/Fs; % observation time
time_vec = (0:Ns-1)*Ts; %Time Vector
freq_vec = (0:Ns-1)*(Fs/Ns);
f_delta = freq_vec(2)-freq_vec(1);

%Create Original Signal
NOF = 201; %Number of Frequencies
if NOF>1
    f0 = (1:NOF)./T0;
    f1 = f0(f0<0.4*Fs);
else
    f1 = 1/T0;
end

amplitudes = sin(30*pi*f1/f1(end))./(30*pi*f1/f1(end));
%phases = 0;
%amplitudes = randn(size(f1));
%phases = 2*pi*rand(size(f1));

original_sig = amplitudes*cos( 2*pi*(f1')*(time_vec-T0/2)) ;

DT = 1e-6;
DD = 100;

%Shifted Signal method 1
int1 = amplitudes*cos( 2*pi*(f1')*(time_vec-DT-T0/2)) ; %Time shift
shifted_sig1=(int1.').*exp(1i*2*pi*DD*time_vec'); %Frequency Shift

%Shifted Signal Method 2
shifted_sig2=time_freq_shift(original_sig,Fs,DT,DD);

%Shifted Signal Method 3
int3=ifft(fftshift((fftshift(fft(original_sig)).*exp(-1i*2*pi*(-Ns/2:Ns/2-1)./Ns*(DT/Ts))))); %Time Shift                
shifted_sig3=int3.'.*exp(1i*2*pi*DD*time_vec'); %Frequency Shift   
    
%Shifted Back
int4 = time_freq_shift(shifted_sig2,Fs,-DT,-DD);

figure;
plot(time_vec*1e6,real(original_sig),'-b',time_vec*1e6, real(shifted_sig1),'-r',...
    time_vec*1e6,real(shifted_sig2),'-g',time_vec*1e6,real(shifted_sig3),'m-',...
    time_vec*1e6,real(int4),'-k');
legend('Original', 'Shifted1', 'Shifted2','Shifted3','Shifted Back');
xlabel 'Time[\mu sec]';
grid
drawnow;

