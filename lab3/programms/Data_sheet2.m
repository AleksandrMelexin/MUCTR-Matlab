%��������� �������� ��������� ������
%������������ �������, �������� �������
%��������� ��������, �������������� ��
%����������� ������ �� ������� 0 �
%����������� ����������� 1
%������� ������� ������������
clear all
%���������� ������������ �������
%������������� ������
nmax=300;
k=0;
%���������� ���� ������ ������������
%������� A - detA
for n=1:5:nmax
    k=k+1;
    order(k)=n;
    %��������� �������� ������� A
    A=randn(n);    
    %��������� ����������� ������� A
    d(k)=det(A);
    %��������� � ��������������� �����
    %��� �������� �������� ������������ 
    d(k)=sign(d(k))*log10(d(k));
end
%������ ������ ����������� ��������
%������������ ������� �� ������� �������
plot(order,d);
norm(A,1)
cond(A,1)
norm(A,2)
cond(A,2)