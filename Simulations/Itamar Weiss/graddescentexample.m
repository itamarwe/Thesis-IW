%An example of gradient descent
%for z=x^2+y^2;
close all;
clear all;

DELTA_X = 2;
DELTA_Y = 2;
MIU = 0.1;

x = rand(1)*100;
y = rand(1)*100;

gradX = inf;
gradY = inf;

histX = [x];
histY = [y];

while (norm([gradX gradY])>0.1)
    gradX = ((x+DELTA_X+randn(1)*0).^2-(x.^2+randn(1)*0))/DELTA_X;
    gradY = ((y+DELTA_Y+randn(1)*0).^2-(y.^2+randn(1)*0))/DELTA_Y;
    x = x-MIU*gradX;
    y = y-MIU*gradY;
    histX = [histX x];
    histY = [histY y];
end

plot3(histX, histY, histX.^2+histY.^2);

