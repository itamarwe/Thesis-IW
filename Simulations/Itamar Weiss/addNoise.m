function [receivedSignal] = addNoise(cleanReceivedSignal,SNRdB,N);
[L,Ntilde] = size(cleanReceivedSignal);
noise = zeros(L,Ntilde);
%Receiver Noise
for l=1:L
    a = sqrt(cleanReceivedSignal(l,:)*cleanReceivedSignal(l,:)');
    b = (10^(-SNRdB/20));
    noise(l,:) = a*b*(1/sqrt(2))*(randn(1,Ntilde)+i*randn(1,Ntilde))/sqrt(Ntilde);
    receivedSignal = cleanReceivedSignal+noise;
end