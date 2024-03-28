p=[2  -8  8  0  -1]; 
roots(p)
%Найдем графическое решение
x=-1:0.1:3;
y=polyval(p,x);
plot(x,y,'-k'),grid