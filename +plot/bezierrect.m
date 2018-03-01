function P = bezierrect(x,y,w,h,c)
% generate points 
xx = [x+c*w,x,x,x+c*w,x+w-c*w,x+w,x+w,x+w-c*w,x+c*w,x];
yy = [y,y+c*h,y+h-c*h,y+h,y+h,y+h-c*h,y+c*h,y,y,y+c*h];
P = spcrv([xx;yy],3)';


end