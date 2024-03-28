%��������� ������ ����� ���������
%������� ������� ��������
%������� ������� ������������
clear all
%���������� �������� ��������� ������ ����� 
eps=1e-10;
%���������� ����������� ��������� �����
%��������
itermax=100;
%���������� �������� ��������� f(x)=x^2-a=0
a=2;
%������ ������� phi
phi=@(x)0.5*(x+a/x);
%���������� ��������� ����������� x0
x0=0.5;
x1=phi(x0); x2=phi(x1);
root(1)=x0; root(2)=x1; root(3)=x2;
it=4;
%���������� ���� ������� ���������������� ��������
while ((x2-x1)^2/abs(2*x1-x0-x2)>eps)&(it<itermax)
    x0=x1; x1=x2; x2=phi(x1);
    root(it)=x2;
    it=it+1;
end
root(it)=phi(x2);
%������ ���������� ���������� �������� �
%����� ���������
plot(1:it,root,'-*');