
function [sig,int, int0, int0_dot, Delayed_sig, Delayed_sig0] = CreateRxSignal(rTransmitter, v, C, Fc, rReceiverMat,B,Fs,Ns,W)
%function [sig,int, int0, int0_dot, Delayed_sig, Delayed_sig0] = CreateRxSignal(rTransmitter, v, C, Fc, rReceiverMat,B,Fs,Ns,W)
%Creates the clean received signal

L = size(rReceiverMat,1);
Ts = 1/Fs;
T0=Ns/Fs; % observation time

time_vec = (0:Ns-1)*Ts; %Time Vector
freq_vec = (0:Ns-1)*(Fs/Ns);
f_delta = freq_vec(2)-freq_vec(1);

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

