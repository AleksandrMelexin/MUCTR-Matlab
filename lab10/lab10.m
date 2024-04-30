clc; clear;
f = @(x) x.^2 + exp(x+3); % �������� �������
f2 = @(x) 2 + exp(x+3); % 2-�� �����������
f4 = @(x) exp(x+3); % 4-�� �����������
% �������� ��������������
a = -5;
b = 5;

tol = 0.01; % �������� ���������� ���������

% ����� ������������� �������� 2-�� �����������
m2 = NaN;
for i=a:0.1:b
    if isnan(m2)
        m2 = f2(i);
    end
    if f2(i) > m2
        m2 = f2(i);
    end    
end

h = sqrt((12 * tol)/(b-a)*m2); % ������ ���� h �� ������� 5
fprintf('��� h ��� ������ ��������: %.10f\n', h);
fprintf('������������ �������� 2-�� ����������� �� ��������� (a; b): %.10f\n', m2);
I_trapezoidal = (h/2) * (f(a) + f(b)); % ����� �������� � ���������� ������� �� ��������� �����
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
plot(h_values, err_values);
xlabel('��� h');
ylabel('�����������');
title('����������� �������� �������������� ������� �������� �� ���� h');

% ����� ������������� �������� 2-�� �����������
m4 = NaN;
for i=a:0.1:b
    if isnan(m4)
        m4 = f4(i);
    end
    if f4(i) > m4
        m4 = f4(i);
    end    
end

h = sqrt((180 * tol)/(b-a)*m4); % ������ ���� h �� ������� 5
fprintf('��� h ��� ������ ��������: %.10f\n', h);
fprintf('������������ �������� 4-�� ����������� �� ��������� (a; b): %.10f\n', m4);
I_simpson = h/3 * (f(a) + 4*f((a+b)/2) + f(b)); % ����� ��������
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

I_matlab = integral(f, a, b); % ������������� ����������� ������� MATLAB ��� ���������� ���������

% ����� �����������
fprintf('�������� ������� ��������: %.10f\n', I_trapezoidal);
fprintf('�������� ������� ��������: %.10f\n', I_simpson);
fprintf('��������, ����������� � ������� ����������� ������� MATLAB: %.10f\n', I_matlab);

%���������� �������������� � �������������� ����������
syms x a p;
f_2 = (a^x) * exp(-x);
f_3 = (1+x)/((x+a).^(p+1));
res = int(f_2, x);
fprintf('������������ ��������: %s\n', res);
res = int(f_3, x, 0, Inf);
fprintf('������������� ��������: %s\n', res);
