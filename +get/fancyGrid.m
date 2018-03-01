function [x,y,w,h] = fancyGrid(A,ig,og,typ,W)
% generate a fancy relative x y grid spacing
% A = a vector were each element is the number of rectangles in each row
% ni = A[i] where i=[1:end] number of columns in each row
% m = size(A,1)
% ig = size of inner gap the spacing between each rectangle (precentage)
% og = outer gap the offset the grid has (precentage)
% ig = 0.05
% og = 0.05
% A = [3 3 2]
%{
[x,y,w,h] = fancyGrid([1,4],0.05,0.05,'weighted',[0.9,0.1])
    A = 1D vector with number of elments at each row;
    N = sum(A);
    
    returns two N long 1D vectors pertainign the position of each object
    use like this :
    A = [0,8,4,2,4,5];
    [x,y] = fancyGrid(A)
    N = sum(A);c=0;n = max(A);m = length(A);w = 1/n;h = 1/m;
    figure;
    for ii=1:length(A)
       for jj=1:A(ii)
           c=c+1;
           p.(['panel',num2str(ii)]) = axes('position',[x(c),y(c),w,h]);
       end
    end
[x,y,w,h] = fancyGrid([5;5;3;3],0.05,0.05,'weighted',[.02;.35;.02;.5])


%}
m = size(A(:),1);
n = max(A);
N = sum(A);
x = zeros(m,n);
y = zeros(m,n);
wt= zeros(m,1);
ht= zeros(m,1);
w = (1-og*2-ig*(n-1))/n;
h = (1-og*2-ig*(m-1))/m;
shift = zeros(m,1);
switch typ
    case 'same'
        for ii=1:m
            if A(ii)==n; shift(ii) = 0;
            else shift(ii) =(1-A(ii)/n)/2;end
        end
        for ii=1:m
            for jj=1:A(ii)
                x(ii,jj) = shift(ii)+(A(ii)-jj)*w+og+ig*(A(ii)-jj);
                y(ii,jj) = (m-ii)*h+og+ig*(m-ii);
            end
        end
    case 'stretch'
        c = 0;
        for ii=1:m
            wt(ii) = (1-og*2-ig*(A(ii)-1))/A(ii);
            ht(ii) = h;
            for jj=1:A(ii)
                c = c+1;
                x(c) = shift(ii)+(A(ii)-jj)*wt(ii)+og+ig*(A(ii)-jj);
                y(c) = (ii-1)*ht(ii)+og+ig*(ii-1);
            end
        end
        
    case 'weighted'
        W = (W./sum(W))*m;
        for ii=1:m
            wt(ii) = (1-og*2-ig*(A(ii)-1))/A(ii);
            ht(ii) = h*W(ii);
            for jj=1:A(ii)
                x(ii,jj) = (A(ii)-jj)*wt(ii)+og+ig*(A(ii)-jj);
                y(ii,jj) = sum(ht(1:(ii-1)))+og+ig*(ii-1);
            end
        end
        w = repelem(wt',A);
        h = repelem(ht',A);
        x = x';
        y = y';
end

ix = x>0;
x = x(ix)';
y = y(ix)';
end