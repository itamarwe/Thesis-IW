%Plots the TOA scenario figure


close all;
clear all;

c = 300;
a = 220;
b = sqrt(c^2-a^2);
x = a:1500;
y = sqrt(b^2*(x.^2/a^2-1));

%Two receivers at (-300,0) and (300,0) try to locate a target.
%The differnce in distances is 2*a = 200[m]

figure;
hold on;
plot([-300,300],[0,0],'kh');
plot(582.45,500,'kp');
legend('Receivers','Transmitter');
plot(x,y,'--');
plot(x,-y,'--');
plot(-x,y,'--');
plot(-x,-y,'--');
grid on;
xlabel 'x[m]';
ylabel 'y[m]';
title 'TDOA scenario';

