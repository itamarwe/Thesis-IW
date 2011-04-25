function CRBMat = CRB(Fc,c,T,B, SNRdB,sNorm2, sDotNorm2, imCorrssDot, rTransmitter,v, rReceiverMat)
%Known Signals CRB
L = size(rReceiverMat,1);
rTransmitterMat = ones(L,2)*[rTransmitter(1) 0;0 rTransmitter(2)];
rDiffMat = rReceiverMat- rTransmitterMat;
d = sqrt(rDiffMat(:,1).^2+rDiffMat(:,2).^2);
vNorm2 = v*v';

sNorm2 = ones(1,L)*sNorm2;
sDotNorm2 = ones(1,L)*sDotNorm2;
imCorrssDot = ones(1,L)*imCorrssDot;

cosTheta = rDiffMat(:,1)./d;
sinTheta = rDiffMat(:,2)./d;
cosPhi = rDiffMat*v'./(d*sqrt(vNorm2));

SNR = 10^(SNRdB/10)*ones(1,L);
FIM = zeros(4); %(1-x,2-y,3-vx,4-vy];

%Calculation of the Fischer Information Matrix

%[J]x,x
for l=1:L
    FIM(1,1) = FIM(1,1)+SNR(l)*cosTheta(l)^2*((pi*vNorm2/d(l)*cosPhi(l))^2+ ...
        B^2/Fc^2*sDotNorm2(l)/sNorm2(l)+2*pi*B/Fc*vNorm2/d(l)*cosPhi(l)*imCorrssDot(l)/sNorm2(l));
end
FIM(1,1) = 4*T/B*(Fc/c)^2*FIM(1,1);

%[J]y,y
for l=1:L
    FIM(2,2) = FIM(2,2)+SNR(l)*sinTheta(l)^2*((pi*vNorm2/d(l)*cosPhi(l))^2+ ...
        B^2/Fc^2*sDotNorm2(l)/sNorm2(l)+2*pi*B/Fc*vNorm2/d(l)*cosPhi(l)*imCorrssDot(l)/sNorm2(l));
end
FIM(2,2) = 4*T/B*(Fc/c)^2*FIM(2,2);

%[J]vx,vx
FIM(3,3) = 4*(pi*Fc/c)^2*T/B*SNR*(cosTheta.^2);

%[J]vy,vy
FIM(4,4) = 4*(pi*Fc/c)^2*T/B*SNR*(sinTheta.^2);

%[J]x,y = [J]y,x
for l=1:L
    FIM(1,2) = FIM(1,2)+SNR(l)*cosTheta(l)*sinTheta(l)*((pi*vNorm2/d(l)*cosPhi(l))^2+ ...
        B^2/Fc^2*sDotNorm2(l)/sNorm2(l)+2*pi*B/Fc*vNorm2/d(l)*cosPhi(l)*imCorrssDot(l)/sNorm2(l));
end
FIM(1,2) = 4*T/B*(Fc/c)^2*FIM(1,2);
FIM(2,1) = FIM(1,2);

%[J]vx,vy = [J]vy,vx
FIM(3,4) = 4*(pi*Fc/c)^2*T/B*SNR*(cosTheta.*sinTheta);
FIM(4,3) = FIM(3,4);

%[J]x,vx = ([J]vx,x)
for l=1:L
    FIM(1,3)= FIM(1,3)+SNR(l)*cosTheta(l)^2*(pi*(vNorm2/d(l)*cosPhi(l)+1)-B/Fc*imCorrssDot(l)/sNorm2(l));
end
FIM(1,3) = FIM(1,3)*4*pi*T/B*(Fc/c)^2;
FIM(3,1) = FIM(1,3);

%[J]y,vy = ([J]y,vy)
for l=1:L
    FIM(2,4)= FIM(2,4)+SNR(l)*sinTheta(l)^2*(pi*(vNorm2/d(l)*cosPhi(l)+1)-B/Fc*imCorrssDot(l)/sNorm2(l));
end
FIM(2,4) = FIM(2,4)*4*pi*T/B*(Fc/c)^2;
FIM(4,2) = FIM(2,4);

%[J]x,vy = [J]y,vx = [J]vy,x = [J]vx,y
for l=1:L
    FIM(1,4) = FIM(1,4)+SNR(l)*sinTheta(l)*cosTheta(l)*(pi*(vNorm2/d(l)*cosPhi(l)+1)-B/Fc*imCorrssDot(l)/sNorm2(l));
end
FIM(1,4) = FIM(1,4)*4*pi*T/B*(Fc/c)^2;
FIM(4,1)=FIM(1,4);
FIM(2,3) = FIM(1,4);
FIM(3,2) = FIM(1,4);

CRBMat = inv(FIM);



