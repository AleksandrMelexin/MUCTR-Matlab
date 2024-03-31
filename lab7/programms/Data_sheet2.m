%������ ��� ������ ����� ��������������� �������  
function norm
%���������� �������� ������ �����
a=-0.1; b=1.1;
%������������ ���������� ����� ������������ 
nmax=16;
%���� �������� ���������������� ���������
%� ������ ����� ������������ �� 3 �� nmax
for i=3:nmax
    h=1.0/(i-1);
    x=0:h:1;
    nor(i-2)=interpol(x,a,b);
end
%������������ ������� ����������� �����
%��������������� ������� �� ����� �����
%������������
n=3:nmax;
plot(n,nor);
%�������������� ������������
function norm=interpol(x,a,b)
%��������� �������� ��������������� �������
%������ ��� �������� ���������� ����������,
%��������������� �� ����������� ������
y=[];
for i=1:length(x)
    y=[y;randn];
end
%������ ������� ��������� (2) � �������
%������-������� ������������� ��������
coef=vander(x)\y;
%������ ��������������� �������
xv=a:0.01:b; 
phi=polyval(coef,xv);
%���������� �������� ����� ���������������
%�������
norm=max(abs(phi));