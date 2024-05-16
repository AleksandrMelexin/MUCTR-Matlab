clc; clear;
f = @(x) x.^2 + exp(x+3); % �������� �������
f2 = @(x) 2 + exp(x+3); % 2-�� �����������
f4 = @(x) exp(x+3); % 4-�� �����������

% �������� ��������������
a = 0;
b = 10;
tol = 0.01; % �������� ���������� ���������

I_matlab = integral(f, a, b); % ������������� ����������� ������� MATLAB ��� ���������� ���������

% ����� ������������� �������� 2-�� �����������
m2 = NaN;
for i = a:0.1:b
    if isnan(m2)
        m2 = f2(i);
    end
    if f2(i) > m2
        m2 = f2(i);
    end
end

h_trap = sqrt((12 * tol) / ((b - a) * m2)); % ������ ���� h ��� ������ ��������
fprintf('��� h ��� ������ ��������: %.10f\n', h_trap);
fprintf('������������ �������� 2-�� ����������� �� ��������� (a; b): %.10f\n', m2);

% ����� �������� ��� ���������
I_trap_no_runge = (h_trap / 2) * (f(a) + f(b)); 
n_trap = 1;
err_trap = tol + 1;
while err_trap > tol
    h_trap = h_trap / 2;
    x_trap = a + h_trap : h_trap : b - h_trap;
    I_new_trap_no_runge = (h_trap / 2) * (f(a) + 2 * sum(f(x_trap)) + f(b));
    err_trap = abs(I_new_trap_no_runge - I_trap_no_runge);
    I_trap_no_runge = I_new_trap_no_runge;
    n_trap = n_trap + 1;
end
h = 10; % ��������� �������� ���� ��� ������ �������� ��� ���������
h_values_trap_no_runge = [];
err_values_trap_no_runge = [];
while h >= 10^-4
    h_values_trap_no_runge = [h_values_trap_no_runge, h];
    x_trap_no_runge = a + h : h : b - h;
    I_new_trap_no_runge = (h / 2) * (f(a) + 2 * sum(f(x_trap_no_runge)) + f(b));
    err_values_trap_no_runge = [err_values_trap_no_runge, abs(I_new_trap_no_runge - I_matlab)];
    h = h / 2;
end

% ����� �������� � ���������� ������� �� ��������� �����
I_trapezoidal = (h_trap / 2) * (f(a) + f(b)); 
n_trap = 1;
err_trap = tol + 1;
while err_trap > tol
    h_trap = h_trap / 2;
    x_trap = a + h_trap : h_trap : b - h_trap;
    I_new_trap = (h_trap / 2) * (f(a) + 2 * sum(f(x_trap)) + f(b));
    err_trap = abs(I_new_trap - I_trapezoidal);
    I_trapezoidal = I_new_trap;
    n_trap = n_trap + 1;
end
h_values_trap = [];
err_values_trap = [];
h = 10; % ��������� �������� ����
while h >= 10^-4
    h_values_trap = [h_values_trap, h];
    x_trap = a + h : h : b - h;
    I_new_trap = (h / 2) * (f(a) + 2 * sum(f(x_trap)) + f(b));
    err_values_trap = [err_values_trap, abs(I_new_trap - I_matlab)];
    h = h / 2;
end

% ����� ������������� �������� 4-�� �����������
m4 = NaN;
for i = a:0.1:b
    if isnan(m4)
        m4 = f4(i);
    end
    if f4(i) > m4
        m4 = f4(i);
    end
end

h_simp = sqrt((180 * tol) / ((b - a) * m4)); % ������ ���� h ��� ������ ��������
fprintf('��� h ��� ������ ��������: %.10f\n', h_simp);
fprintf('������������ �������� 4-�� ����������� �� ��������� (a; b): %.10f\n', m4);

I_simpson = (h_simp / 3) * (f(a) + 4 * f((a + b) / 2) + f(b)); % ����� ��������
n_simp = 1;
err_simp = tol + 1;

% ���������� ��������� ������� ��������
while err_simp > tol
    h_simp = h_simp / 2;
    x_simp = a + h_simp : 2 * h_simp : b - h_simp;
    I_new_simp = (h_simp / 3) * (f(a) + 4 * sum(f(x_simp)) + 2 * sum(f(x_simp + h_simp)) + f(b));
    err_simp = abs(I_new_simp - I_simpson);
    I_simpson = I_new_simp;
    n_simp = n_simp + 1;
end

% ������ ����������� �������� �� ���� h ��� ������ ��������
h_values_simp = [];
err_values_simp = [];
h = 10; % ��������� �������� ����
while h >= 10^-4
    h_values_simp = [h_values_simp, h];
    x_simp = a + h : 2 * h : b - h;
    I_new_simp = (h / 3) * (f(a) + 4 * sum(f(x_simp)) + 2 * sum(f(x_simp + h / 2)) + f(b));
    err_values_simp = [err_values_simp, abs(I_new_simp - I_matlab)];
    h = h / 2;
end

figure;
loglog(h_values_trap, err_values_trap, '-o', 'DisplayName', '����� �������� � ����������');
hold on;
loglog(h_values_simp, err_values_simp, '-', 'DisplayName', '����� ��������');
loglog(h_values_trap_no_runge, err_values_trap_no_runge, '-x', 'DisplayName', '����� �������� ��� ���������');
xlabel('��� h');
ylabel('�����������');
title('����������� �������� �������������� �� ���� h');
legend show;

% ����� �����������
fprintf('�������� ������� �������� ��� ���������: %.10f\n', I_trap_no_runge);
fprintf('�������� ������� �������� � ����������: %.10f\n', I_trapezoidal);
fprintf('�������� ������� ��������: %.10f\n', I_simpson);
fprintf('��������, ����������� � ������� ����������� ������� MATLAB: %.10f\n', I_matlab);

% ���������� �������������� � �������������� ����������
syms x a p;
f_2 = (a^x) * exp(-x);
f_3 = (1 + x) / ((x + a) .^ (p + 1));
res = int(f_2, x);
fprintf('������������� ��������: %s\n', char(res));
res = int(f_3, x, 0, Inf);
fprintf('������������� ��������: %s\n', char(res));