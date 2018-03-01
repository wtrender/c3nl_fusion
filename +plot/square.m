function [xp,yp] = square(x,y,w,h)
%x and y are the coordinates of the center of the square
% w is the width and h is the height
%
%
x1 = x-w/2;x2=x+w/2;y1 =y-h/2; y2 =y+h/2;
xp= [ x1 x1 x2 x2 x1] ;
yp= [ y1 y2 y2 y1 y1];
end