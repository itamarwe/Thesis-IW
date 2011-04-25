%Simulation 01 - Physical Intuition for 1 known Fc and Velocity with 2
%receivers

clear all;
close all;

C = 3e08; %[m/s]
r1 = [-10000 0]; %Position of the 1st receiver [m]
r2 = [100000 0]; %Postion of the 2nd receiver [m]

r = [0 50000]; %Target position [m]
v = [100/sqrt(2) 100/sqrt(2)];%Target velocity [m/s]
Fc = 900e06; %[Hz]

%True values
psi = atan(v(2)/v(1));
dr1 = r-r1;
dr2 = r-r2;
theta1 = atan(dr1(2)/dr1(1));
theta2 = atan(dr2(2)/dr2(1));

deltaF1 = -Fc*((1/C)*(v*dr1'))/sqrt(dr1*dr1');
deltaF2 = -Fc*((1/C)*(v*dr2'))/sqrt(dr2*dr2');

% Measured Values
mFc = Fc;
mDeltaF1 = -mFc*((1/C)*(v*dr1'))/sqrt(dr1*dr1');
mDeltaF2 = -mFc*((1/C)*(v*dr2'))/sqrt(dr2*dr2');

mPhi1 = acos(-C/sqrt(v*v')*mDeltaF1/mFc);
mTheta1 = psi-mPhi1;

mPhi2 = -acos(-C/sqrt(v*v')*mDeltaF2/mFc);
mTheta2 = psi-mPhi2;

%plot the scene
figure;
hold on;
grid on;
title({'Diffrential Doppler Geolocation','Single Frequency, Known Velocity, Known Carrier'});
plot(r1(1),r1(2),'x');
plot(r2(1),r2(2),'x');
plot(r(1),r(2),'or');
plot([r1(1) r(1)], [r1(2), r1(2) + tan(theta1)*dr1(1)],'b');
plot([r2(1) r(1)], [r2(2), r2(2) + tan(theta2)*dr2(1)],'b');

plot([r1(1) r(1)], [r1(2), r1(2) + tan(mTheta1)*dr1(1)],'r-.');
plot([r2(1) r(1)], [r2(2), r2(2) + tan(mTheta2)*dr2(1)],'r-.');


%if carrier frequency is unknown

mF1 = 1e6+1e-7;%mFc+mDeltaF1;%
mF2 = 1e6-1.5e-7;%mFc + mDeltaF2; 

maxFc = 1/(1-sqrt(v*v')/C)*min([mF1, mF2]);
minFc = 1/(1+sqrt(v*v')/C)*max([mF1,mF2]);

Fc = linspace(minFc,maxFc,1000);

mFc = Fc;

mDeltaF1 = mF1 - mFc;
mDeltaF2 = mF2 - mFc;

mPhi1 = acos((-C/sqrt(v*v')*mDeltaF1)./mFc);
mTheta1 = psi-mPhi1;

mPhi2 = -acos(-C/sqrt(v*v')*mDeltaF2./mFc);
mTheta2 = psi-mPhi2;

mR(1,:) = (r1(1)-r2(1) - (tan(mTheta1)*r1(1)-tan(mTheta2)*r2(1)))./(tan(mTheta1)-tan(mTheta2));

mR(2,:) = r1(2)+(mR(1,:)-r1(1)).*tan(mTheta1);

plot(mR(1,:),mR(2,:));
