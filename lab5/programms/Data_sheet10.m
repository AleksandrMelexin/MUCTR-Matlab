%��������� ������������ ������ ������� fzero
function testfzero
%������� ������� ������������
clear all
global n
%���������� ����� ��������: 1,2,...,n
n=160;
%���������� ������ ��������� �������� ���
%������������� ��������, ������������� ���
%������ ������� fzero
x0=0.5:0.1:n;
%���������� ���� ������������� ������� fzero
for k=1:length(x0)
    y(k)=fzero(@f,x0(k));
end
%������ ������ ����������� ������ ��������,
%������� ���������� ������� fzero �� ���������
%�������� x0 (������ ���� ��������� ��������)
plot(x0,y);
%���������� �������, ����� ������� ������� fzero
function y=f(x)
global n
y=1;
for i=1:n
    y=y*(x-i);
end