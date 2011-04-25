
function [originalSignal, cleanReceivedSignal, sNorm2, sDotNorm2,imCorrssDot] = CreateRxSignal(rTransmitter, v, C, Fc, rReceiverMat,B,Fs,N)
%******************************************************
%[receivedSignal, transmittedSignal] = CREATE_RX_SIGNAL(B, Fs,
%N(+maxGridSampleDelay),rTransmitter,v,rReceiverMat)
%****************************************************8
STEPSIGNAL = 1;

L = size(rReceiverMat,1);
rTransmitterMat = ones(L,2)*[rTransmitter(1) 0;0 rTransmitter(2)];
rDiffMat = rTransmitterMat-rReceiverMat;
rDistances = sqrt(rDiffMat(:,1).^2+rDiffMat(:,2).^2);

%Calculate the delay for the true position of the transmitter
%Get the delays in [Sec]
[maxDelay delay] = CalcMicroSecRxSignalDelay(rTransmitter, rReceiverMat,C);
%[sampleMaxDelay sampleDelay]=CalcRxSignalDelay(rTransmitter,%rReceiverMat,C,Fs);

%Original Signal
%s = randn(1,N+sampleMaxDelay).*exp(j*randn(1,N+sampleMaxDelay));
%s = sin(2*pi*Fsig*(1:(N+sampleMaxDelay))/Fs);
%[k,n] = meshgrid(0:(N-1),0:(N-1));
%kn = zeros(N);
%ak = 1/sqrt(2)*(randn(1,N)+i*randn(1,N));

%kn = k.*n;

if (STEPSIGNAL)
    %The width in samples of the pulse(must be an odd number)
W = 2*100+1;
k = 0:N-1;

%Fourier transform of a 1/W high step with width W
ak = sin(pi*W*k/N)./(W*sin(pi*k/N));
ak(1) = 1;

%create a N/2 delay so that the pulse will be in the middle of the interval
m = floor(N/2);
ak = ak.*exp(-2*pi*i/N*k*m);
end

%Band Limiting the signal

Kmax= floor(B/2*N/Fs);
ak(Kmax+2:N-Kmax)=0;

%Normalize ak so |s|^2=1
sNorm = sqrt(ak*ak'/N);
ak = ak/sNorm;

%Find the coefficients for sDot
akDot = ak.*(2*pi*i*Fs/N*k);

%Reconstruction

[kMat,nMat] = meshgrid(0:N-1,0:N-1);
F = exp(2*pi*i/N*kMat.*nMat);
s = 1/N*ak*F;
sDot = 1/N*akDot*F;

originalSignal = s;


sNorm2 = sNorm^2;

imCorrssDot = imag(sDot*s');

sDotNorm2 = akDot*akDot'/N;


%Complex Path Attenuation
b = exp(i*2*pi*rand(1,L));
%b=ones(1,L);

%Delayed Signal
delayedSignal = zeros(L,N);

for j=1:L
    m = delay(j)*(Fs*1e-6);
    bk = ak.*exp(-2*pi*i/N*k*m);
    delayedSignal(j,:) = 1/N*bk*F;
end %j

%A = zeros(N,N,L);

cleanReceivedSignal = zeros(L,N);
miu = zeros(L,1);
hzDopplerShift = zeros(L,1);
%Doppler Shift
for l=1:L
    %miu calculation
    miu(l) = -1/C*v*rDiffMat(l,:)'/rDistances(l);
    %frequency calculation
    hzDopplerShift = Fc*miu(l);
    %Al Matrix
    %A(:,:,l) = diag(exp(2*pi*i*hzDopplerShift(l)*(1:N)/Fs));
    vA = exp(2*pi*i*hzDopplerShift*(1:N)/Fs);
    cleanReceivedSignal(l,:) =b(l)*vA.*delayedSignal(l,:);
end


% Received Signal
% cleanReceivedSignal = zeros(L,N);
% 
% for l=1:L
%     vA=exp(2*pi*i*hzDopplerShift(l)*(1:N)/Fs);
%     cleanReceivedSignal(l,:) = b(l)*vA.*delayedSignal(l,:);
% end




