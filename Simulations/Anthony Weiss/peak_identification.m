function [peak] = peak_identification(x, y, I)
y1 = y(I-1);
y2 = y(I);
y3 = y(I+1);
delta = x(2)-x(1);
if abs(y1+y3-2*y2)>eps*5
    temp = (y1-y3)*delta/2/(y1+y3-2*y2);
else
    temp =0;
end
peak = x(I) + temp;