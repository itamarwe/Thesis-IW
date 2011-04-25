%Plots the TOA scenario figure

close all;
clear all;

%Two receivers at (-300,0) and (300,0) try to locate a target.

angles =0:0.01:2*pi;

figure;
hold on;
plot([-300,300],[0,0],'kh');
plot(364.5833,747.2406,'kp');
legend('Receivers','Transmitter');
plot(1000*cos(angles)-300,1000*sin(angles),'--');
plot(750*cos(angles)+300,750*sin(angles),'--');
grid on;
xlabel 'x[m]';
ylabel 'y[m]';
title 'TOA scenario';

