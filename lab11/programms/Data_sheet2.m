%����������� �.�. �������������� ������.
%������ � �������� � ����� MATLAB
%��������� ���������� �������������
%���������������� �������� � ������ 
%������� ���������� ������������� �����
%������� ������� ������������
clear all
%���������� ������������� ����������
syms xs h
%���������� �������� ��������������� �������
%� ������ xs-h, xs � xs+h ��������������
syms Fsm Fs Fsp
%���������� ������� �������� �������
%��������� (14) ������������ ����������� a, b, c
A=[1 xs-h (xs-h)^2;1 xs xs^2;1 xs+h (xs+h)^2];
%���������� ������ ����� r ��������
%������� ���������
r=[Fsm;Fs;Fsp];
%������������ ������ �������� ������� ���������
abc=A\r;
%������� ��������� ����� ���������� �������� xps
%���������, ��� ��������, �� ����������� -b/2c
xps=-abc(2)/(2*abc(3))
