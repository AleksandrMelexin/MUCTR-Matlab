%������������ � ������� ����������� ������� �������
%��������� ������ ����� ������������
x=0:0.025:1;
%��������� �������� ��������������� �������,
%������ ��� �������� ���������� ����������,
%��������������� �� ����������� ������
y=[];
for i=1:length(x)
    y=[y randn];
end
%������ ����� �� ������� ������������
xv=0:0.001:1.0;
%���������� � ����������� ��������� MATLAB
yv=interp1(x,y,xv,'cubic');
%������ ������
plot(x,y,'*',xv,yv);