[z,f,exitflag,output]=fminsearch(@(x) sqrt(x(1)^2+x(2)^2), [2,2])
%Построение графика
[x y]=meshgrid(-2:0.2:2, -2:0.2:2);
z=sqrt(x.^2+y.^2);
surf(x,y,z);