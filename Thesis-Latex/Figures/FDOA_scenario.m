%Plots the TOA scenario figure


close all;
clear all;

p = [0;300];

p1 = [-300;0];
p2 = [300;0];
v1 = [10;0];
v2 = [10;0];


%The difference between the relative radial velocity: deltaF12*c/Fc
deltaVr12 = (p-p1)'*v1/norm(p-p1)-(p-p2)'*v2/norm(p-p2);

theta2min = acos(min(1,(norm(v1)-deltaVr12)/norm(v2)));
theta2max = acos(max(-1,(-norm(v1)-deltaVr12)/norm(v2)));

theta2 = linspace(theta2min,theta2max,100);

%theta2 = atan(v2(2)/v2(1))-atan((p(2)-p2(2))/(p(1)-p2(1)));
%theta2 = acos(v2'*(p-p2)/(norm(p-p2)*norm(v2)));

%theta1_known = atan(v1(2)/v1(1))-atan((p(2)-p1(2))/(p(1)-p1(1)));

%theta2 = theta2+linspace(-pi,pi/2,100);

%a unit vector in the direction of p2-p1
ihat = (p2-p1)/norm(p2-p1);
%a unit vector perpendicular to the direction of p2-p1
jhat = [-ihat(2);ihat(1)];

dif = p2-p1;
%The angle between dif and the x-axis
difAngle = acos(dif(1)/norm(dif))*sign(dif(2));
%The distance between the receivers
B = norm(dif);

theta1 = acos((deltaVr12+norm(v2)*cos(theta2))/norm(v1));
theta1 = [theta1 fliplr(-theta1)];
theta2 = [theta2 fliplr(-theta2)];
%theta1 = [theta1 fliplr(-theta1)];
%theta2 = [theta2 fliplr(pi-theta2)];
v1Angle = acos(v1(1)/norm(v1))*(sign(v1(2)));
v2Angle = acos(v2(1)/norm(v2))*(sign(v2(2)));
alpha1 = v1Angle-theta1-difAngle;
alpha2 = v2Angle-theta2-difAngle;




ri = B*tan(alpha1)./(tan(alpha2)-tan(alpha1));
rj = ri.*tan(alpha2);

rx = p2(1)+ihat(1)*ri+jhat(1)*rj;
ry = p2(2)+ihat(2)*ri+jhat(2)*rj;

figure;
hold on;
plot([p1(1),p2(1)],[p1(2),p2(2)],'kh');
plot(p(1),p(2),'kp');
legend('Receivers','Transmitter');
plot(rx,ry,'--');
%plot(rx,-ry,'--');
grid on;
xlabel 'x[m]';
ylabel 'y[m]';
title 'FDOA scenario';

