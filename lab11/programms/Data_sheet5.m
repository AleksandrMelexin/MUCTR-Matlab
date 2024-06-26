%����������� �.�. �������������� ������.
%������ � �������� � ����� MATLAB
%���������, �������������� ����� ��������
%������� F(x,y)=(y-sin(x))^2+0.1x^2
%������� ��������������� ������
%������� ������� ������������
clear all
%������ �������� �������� ������� �����������
%������� F � ����
eps=1e-3;
%���������� ������� � �� ������� �����������
F=@(x,y)(y-sin(x))^2+0.1*x^2;
Fx=@(x,y)-2*cos(x)*(y-sin(x))+0.2*x;
Fxx=@(x,y)2*y*sin(x)+2*cos(2*x)+0.2;
Fy=@(x,y)2*(y-sin(x));
Fyy=@(x,y)2;
%������ ��������� �����������
x(1)=4.5; y(1)=-1;
%������ ������� ����� ����� � ������ ������ �
%������������ ����� �����
k=1; iterm=200;
%���������� ���� ��������������� ������
while ((abs(Fx(x(k),y(k)))>eps)|...
       (abs(Fy(x(k),y(k)))>eps))&(k<iterm)
    %���� ������ �� ���������� x ���������� �
    %������� ������ �������
    xs=x(k);
    for i=1:2
        xs=xs-Fx(xs,y(k))/Fxx(xs,y(k));
    end
    k=k+1;
    x(k)=xs; y(k)=y(k-1);
    %���� ������ �� ���������� y ���������� �
    %������� ������ �������
    ys=y(k);
    for i=1:2
        ys=ys-Fy(x(k),ys)/Fyy(x(k),ys);
    end
    k=k+1;
    x(k)=x(k-1); y(k)=ys;
end
%���������������� ����������� � ����������
%����� ������
[u v]=meshgrid(-2*pi:0.1:2*pi,-pi:0.1:pi);
Func=(v-sin(u)).^2+0.1*u.^2;
%���������� �������� �������, ����� ������
%������� ����� ���������
for i=1:5
    s(i)=F(x(i),y(i));
end
%���������� ����� ������
contour(u,v,Func,s);
hold on
%���������� ���������� ������ � ��������
%������� F
line(x,y,'Color','black');