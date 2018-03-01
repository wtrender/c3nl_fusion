function [x,y]= arc(r1,r2,pc)

% arc on circle
% n = 100; % The number of points in the arc
% v1 = p1-pc;
% v2 = p2-pc;
% c = det([v1,v2]); % "cross product" of v1 and v2
% a = linspace(0,atan2(abs(c),dot(v1,v2)),n); % Angle range
% v3 = [0,-c;c,0]*v1; % v3 lies in plane of v1 and v2 and is orthog. to v1
% v = v1*cos(a)+((norm(v1)/norm(v3))*v3)*sin(a); % Arc, center at (0,0)
% x = v(1,:)+pc(1);
% y= v(2,:)+pc(2);
% 
%  pts = f(t,p1,pc,p2)
% shallow clockwise arc
%     a = atan2(-(x2-x1),y2-y1); % Perpendicular bisector angle
%     b = asin(d/2/r); % Half arc angle
%     c = linspace(a-b,a+b); % Arc angle range
%     e = sqrt(r^2-d^2/4); % Distance, center to midpoint
%     x = (x1+x2)/2-e*cos(a)+r*cos(c); % Cartesian coords. of arc
%     y = (y1+y2)/2-e*sin(a)+r*sin(c);
if ~exist('pc','var');pc = [0,0];end % assume center is in 0,0
x1 = cos(r1);x2 = cos(r2);
y1 = sin(r1);y2 = sin(r2);
t = linspace(0,1,101);
r =  sqrt((x2-pc(1))^2+(y2-pc(2))^2);
d = sqrt((x2-x1)^2+(y2-y1)^2);% Distance between points
f = @(t,xy1,cxy,xy2) kron((1-t).^2,xy1) + kron(2*(1-t).*t,cxy) +  kron(t.^2,xy2); % define a second order Bernstein polynomial

if r/d>1

    
    u  = [cos(r1);sin(r1)];
    v  = [cos(r2);sin(r2)];
    x0 = -(u(2)-v(2))/(u(1)*v(2)-u(2)*v(1));
    y0 =  (u(1)-v(1))/(u(1)*v(2)-u(2)*v(1));
    r  = sqrt(x0^2 + y0^2 - 1);
    thetaLim(1) = atan2(u(2)-y0,u(1)-x0);
    thetaLim(2) = atan2(v(2)-y0,v(1)-x0);
    if u(1) >= 0 && v(1) >= 0
        % ensure the arc is within the unit disk
        theta = [linspace(max(thetaLim),pi,50),...
            linspace(-pi,min(thetaLim),50)].';
    else
        theta = linspace(thetaLim(1),thetaLim(2)).';
    end
    x = r*cos(theta)+x0;
    y = r*sin(theta)+y0;


    
else 
    
    [x,y] = c3nl.deal(f(t',[x1,y1],pc,[x2,y2]));
    
end

end