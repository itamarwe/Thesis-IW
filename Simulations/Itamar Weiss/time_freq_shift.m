function shifted_sig = time_freq_shift(sig,Fs,DT,DD)
%shifted_sig = time_freq_shift(sig,Fs, DT,DD);
%A function that shifts the signal by DT in time and DD in frequency

Ns = length(sig);
Ts = 1/Fs;
time_vec = (0:Ns-1)*Ts; %Time Vector

 if size(sig,1)<size(sig,2)
     int1 = sig.';
 else
     int1=sig;
 end


if (DD~=0)
    int1=int1.*exp(1i*2*pi*DD*time_vec'); %Frequency Shift
end

if (DT~=0)
shifted_sig=ifft(fftshift((fftshift(fft(int1.')).*exp(-1i*2*pi*(-Ns/2:Ns/2-1)./Ns*(DT/Ts))))); %Time Shift                
else
    shifted_sig=int1;
end

 if size(sig,1)<size(sig,2)
     shifted_sig = shifted_sig.';
 end
