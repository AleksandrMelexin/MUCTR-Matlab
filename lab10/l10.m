clc; clear;

f = @(x) cos(x) + x.^2;
f2 = @(x) -cos(x) + 2;
f4 = @(x) cos(x);
% �������� ��������������
a = 0;
b = 1;

% �������� ���������� ���������
tol = 1e-2;

% ����� �������� � ���������� ������� �� ��������� �����
%h = (b - a);
m2 = NaN;
for i=a:0.1:b
    if isnan(m2)
        m2 = f2(i);
    end
    if f2(i) > m2
        m2 = f2(i);
    end    
end   


h = sqrt((12 * tol)/(b-a)*m2);
I_trapezoidal = (h/2) * (f(a) + f(b));
n = 1;
err = tol + 1;
while err > tol
    h = h/2;
    x = a + h : h : b - h;
    I_new = (h/2) * (f(a) + 2 * sum(f(x)) + f(b));
    err = abs(I_new - I_trapezoidal);
    I_trapezoidal = I_new;
    n = n + 1;
end


% ������ ����������� �������� �� ���� h
h_values = zeros(1, n);
err_values = zeros(1, n);
h = (b - a);
for i = 1 : n
    h_values(i) = h;
    x = a + h : h : b - h;
    I_new = (h/2) * (f(a) + 2 * sum(f(x)) + f(b));
    err_values(i) = abs(I_new - I_trapezoidal);
    h = h/2;
end

figure;
plot(h_values, err_values, '-o');
xlabel('��� h');
ylabel('�����������');
title('����������� �������� �������������� ������� �������� �� ���� h');

% ����� ��������

tol = 1e-4;

m4 = NaN;

for i=a:0.1:b
    if isnan(m4)
        m4 = f4(i);
    end
    if f4(i) > m4
        m4 = f4(i);
    end    
end

%h = (b - a) / 2;
h = sqrt((180 * tol)/(b-a)*m4);
I_simpson = h/3 * (f(a) + 4*f((a+b)/2) + f(b));
n = 1;
err = tol + 1;
while err > tol
h = h/2;
x = a + h : 2* h : b - h;
I_new = h/3 * (f(a) + 4 * sum(f(x)) + 2 * sum(f(x + h)) + f(b));
err = abs(I_new - I_simpson);
I_simpson = I_new;
n = n + 1;
end


% ������������� ���������� ��������� ����������� �������
%f4 = @(x) 12 * exp(-x.^2) * (4 * x.^4 - 12 * x.^2 + 3);
%max_f4 = max(f4(a:b));

% ������������� ����������� ������� MATLAB ��� ���������� ���������
I_matlab = integral(f, a, b);

% ����� �����������
fprintf('�������� ������� ��������: %.10f\n', I_trapezoidal);
fprintf('�������� ������� ��������: %.10f\n', I_simpson);
%fprintf('������������ �������� 4-� ����������� �������: %.10f\n', max_f4);
fprintf('��������, ����������� � ������� ����������� ������� MATLAB: %.10f\n', I_matlab);