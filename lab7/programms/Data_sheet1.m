%����������� �.�. �������������� ������.
%������ � �������� � ����� MATLAB
%������������ �������������� ������������
%��������� ������ ����� ������������
x=[0 0.1 0.2 0.3 0.35 0.6 0.7 0.9 0.95 1];
%��������� �������� ��������������� �������
%������ ��� �������� ���������� ����������,
%��������������� �� ����������� ������
y=[];
for i=1:length(x)
    y=[y;randn];
end
%������ ������� ��������� (2) � �������
%������-������� ������������� ��������
a=vander(x)\y;
%������ ��������������� �������
xv=-0.01:0.01:1.01; 
phi=polyval(a,xv);
%������ ������ ��������������� �������,
%� ��������� ����� �������� ������� �
%����� ������������
plot(x,y,'*',xv,phi);