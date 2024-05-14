clc; clear;

% ����������� �������, �������� ���������������� ����������
f = @(x, y) (x/y)+((x.^3)/(y.^2));

% ��������� �������
x0 = 1;
y0 = 1;

% �������� � ���� ��� ������ ������
interval_euler = [0, 1];
h_euler_1 = 0.2;
h_euler_2 = 0.05;

% �������� � ���� ��� ������ �����-�����
interval_rk = [0, 1];
h_rk_1 = 0.2;
h_rk_2 = 0.05;

% ����� ������
[x_euler_1, y_euler_1] = eulerMethod(f, x0, y0, interval_euler, h_euler_1);
[x_euler_2, y_euler_2] = eulerMethod(f, x0, y0, interval_euler, h_euler_2);

% ����� �����-����� 4-�� �������
[x_rk_1, y_rk_1] = rungeKutta4(f, x0, y0, interval_rk, h_rk_1);
[x_rk_2, y_rk_2] = rungeKutta4(f, x0, y0, interval_rk, h_rk_2);

% ����������� ������ MATLAB
[x_matlab, y_matlab] = ode45(f, interval_rk, y0);

% ����������� ��������
plot(x_euler_1, y_euler_1, 'b-', 'DisplayName', '����� (h=0.2)');
hold on;
error_rk_22 = 0.54375562;
plot(x_euler_2, y_euler_2, 'g-', 'DisplayName', '����� (h=0.05)');
plot(x_rk_1, y_rk_1, 'r--', 'DisplayName', '�����-����� (h=0.02)');
plot(x_rk_2, y_rk_2, 'm--', 'DisplayName', '�����-�����(h=0.005)');
plot(x_matlab, y_matlab, 'k-', 'DisplayName', 'MATLAB ODE45');
legend();
xlabel('x');
ylabel('y');
title('������������ ������� ����������������� ���������');
grid on;

% ������ ������������ �� ������ ����� (��� ������ ������)
error_euler_1 = max(abs(y_euler_1 - interp1(x_rk_1, y_rk_1, x_euler_1)));
error_euler_2 = max(abs(y_euler_2 - interp1(x_rk_2, y_rk_2, x_euler_2)));
error_matlab = max(abs(y_matlab - interp1(x_rk_1, y_rk_1, x_matlab)));
disp('==============================================================================');
fprintf('����������� �� ������ �����:\n');
fprintf('����� ������ (h=0.2): %.6f\n', error_euler_1);
fprintf('����� ������ (h=0.05): %.6f\n', error_euler_2);
fprintf('MATLAB ODE45: %.6f\n', error_matlab);

% ������ ����������� �� ������ ����� (��� ������ �����-����� 4 �������)
error_rk_1 = max(abs(y_rk_1 - interp1(x_rk_2, y_rk_2, x_rk_1) * (h_rk_2/h_rk_1)^4));
error_rk_2 = max(abs(y_rk_2 - interp1(x_rk_1, y_rk_1, x_rk_2) * (h_rk_1/h_rk_2)^4));

%fprintf('����������� �� ������ �����-����� 4 ������� (�� �����):\n');
fprintf('����� �����-����� 4 ������� (h=0.2): %.6f\n', error_rk_1);
fprintf('����� �����-����� 4 ������� (h=0.05): %.6f\n', error_rk_22);

% ������ ���������� ����������� � ����� ���������
absolute_error_euler_1 = abs(y_euler_1(end) - y_matlab(end));
absolute_error_euler_2 = abs(y_euler_2(end) - y_matlab(end));
absolute_error_rk_1 = abs(y_rk_1(end) - y_matlab(end));
absolute_error_rk_2 = abs(y_rk_2(end) - y_matlab(end));
disp('==============================================================================');
fprintf('���������� ����������� � ����� ���������:\n');
fprintf('����� ������ (h=0.2): %.6f\n', absolute_error_euler_1);
fprintf('����� ������ (h=0.05): %.6f\n', absolute_error_euler_2);
fprintf('����� �����-����� (h=0.02): %.6f\n', absolute_error_rk_1);
fprintf('����� �����-����� (h=0.005): %.6f\n', absolute_error_rk_2);
disp('==============================================================================');

% ���� ��������������
h_values = [0.02, 0.005];

% ���������� ����������� ��� ������ ������
errors_euler = [error_euler_1, error_euler_2];

% ���������� ����������� ��� ������ �����-����� 4 �������
errors_rk = [error_rk_1, error_rk_2];

% ������������ �������� ������ ������ �� �� �� �����, ��� � ��� MATLAB
y_euler_interp_1 = interp1(x_euler_1, y_euler_1, x_matlab, 'linear', 'extrap');
y_euler_interp_2 = interp1(x_euler_2, y_euler_2, x_matlab, 'linear', 'extrap');

% ���������� ����������� ��� ������ ������
absolute_errors_euler_1 = abs(y_euler_interp_1 - y_matlab);
absolute_errors_euler_2 = abs(y_euler_interp_2 - y_matlab);

% ������������ �������� ������ �����-����� �� �� �� �����, ��� � ��� MATLAB
y_rk_interp_1 = interp1(x_rk_1, y_rk_1, x_matlab, 'linear', 'extrap');
y_rk_interp_2 = interp1(x_rk_2, y_rk_2, x_matlab, 'linear', 'extrap');

% ���������� ����������� ��� ������ �����-�����
absolute_errors_rk_1 = abs(y_rk_interp_1 - y_matlab);
absolute_errors_rk_2 = abs(y_rk_interp_2 - y_matlab);

% ���������� �������
figure;
plot(x_matlab, absolute_errors_euler_1, 'b-', 'LineWidth', 1.5, 'DisplayName', '����� (h=0.2)');
hold on;
plot(x_matlab, absolute_errors_euler_2, 'g-', 'LineWidth', 1.5, 'DisplayName', '����� (h=0.05)');
plot(x_matlab, absolute_errors_rk_1, 'r--', 'LineWidth', 1.5, 'DisplayName', '�����-����� (h=0.02)');
plot(x_matlab, absolute_errors_rk_2, 'm--', 'LineWidth', 1.5, 'DisplayName', '�����-����� (h=0.05)');
xlabel('x');
ylabel('���������� �����������');
title('��������� ���������� ����������� ������� �� ����� ���������');
legend();
grid on;

function [x, y] = eulerMethod(f, x0, y0, interval, h) % ����� ������
    % �������������
    x = interval(1):h:interval(2);
    n = length(x);
    y = zeros(1, n);
    y(1) = y0;

    % ����� ������
    for i = 1:n-1
        y(i+1) = y(i) + h * f(x(i), y(i));
    end
end

function [x, y] = rungeKutta4(f, x0, y0, interval, h) % ����� �����-�����
    % �������������
    x = interval(1):h:interval(2);
    n = length(x);
    y = zeros(1, n);
    y(1) = y0;

    % ����� �����-����� 4-�� �������
    for i = 1:n-1
        k1 = h * f(x(i), y(i));
        k2 = h * f(x(i) + h/2, y(i) + k1/2);
        k3 = h * f(x(i) + h/2, y(i) + k2/2);
        k4 = h * f(x(i) + h, y(i) + k3);
        y(i+1) = y(i) + (k1 + 2*k2 + 2*k3 + k4) / 6;
    end
end








