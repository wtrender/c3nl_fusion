function [X,Y,Z] = tube(xyz,S,r,s)
sz = size(xyz);
if sz(1)~=3;disp('this function is intended to work with 3D points where each row corrosponds to a dimension');return;end
if sz(2)<2;disp('this function needs at least 2 points to work');return;end

if sz(2) > 2 
    n = 100;
    s = 10;
    xyz = interp3dCurve(xyz,n)';
    t = diffShift(xyz);
    n = diffShift(t);
    b = cross(t,n,2);
    b = b./repmat(mnorm(b),1,size(b,2));
    theta=(0:(2*pi/(s-1)):(2*pi));
    X = m_add(xyz(:,1) ,m_mult(r(:),(m_mult(n(:,1),cos(theta)) + m_mult(b(:,1),sin(theta)))));
    Y = m_add(xyz(:,2) ,m_mult(r(:),(m_mult(n(:,2),cos(theta)) + m_mult(b(:,2),sin(theta)))));
    Z = m_add(xyz(:,3) ,m_mult(r(:),(m_mult(n(:,3),cos(theta)) + m_mult(b(:,3),sin(theta)))));
else
    p1 = xyz(:,1);
    p2 = xyz(:,2);
    v=norm(p1-p2);% length of single vector
    t=linspace(0,2*pi,S)'; % span
    Vx=[1 0 0]; % unit vector in x axis
    alpha=acos(dot(Vx,(p1-p2))/(norm(Vx)*v))*180/pi;
    azel=cross(Vx,p1-p2);
    X = repmat([0 v],S,1);
    Y = [r*cos(t),r*cos(t)];
    Z = [r*sin(t),r*sin(t)];
    u = azel(:)/norm(azel);% direction vector
    alph = alpha*pi/180;
    cosa = cos(alph);
    sina = sin(alph);
    vera = 1 - cosa;
    [x,y,z] = c3nl.deal(u');
    rot = [cosa+x^2*vera x*y*vera-z*sina x*z*vera+y*sina; ...
       x*y*vera+z*sina cosa+y^2*vera y*z*vera-x*sina; ...
       x*z*vera-y*sina y*z*vera+x*sina cosa+z^2*vera]';
    [m,n] = size(X);
    newxyz = [X(:), Y(:), Z(:)];
    newxyz = newxyz*rot;   
    X = p2(1)+reshape(newxyz(:,1),m,n);
    Y = p2(2)+reshape(newxyz(:,2),m,n) ;
    Z = p2(3)+reshape(newxyz(:,3),m,n);
end
end
