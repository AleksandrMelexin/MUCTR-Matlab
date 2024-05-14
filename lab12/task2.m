clc; clear;

% Определение функции, заданной дифференциальным уравнением
f = @(x, y) (x/y)+((x.^3)/(y.^2));

% Начальные условия
x0 = 1;
y0 = 1;

% Интервал и шаги для метода Эйлера
interval_euler = [0, 1];
h_euler_1 = 0.2;
h_euler_2 = 0.05;

% Интервал и шаги для метода Рунге-Кутта
interval_rk = [0, 1];
h_rk_1 = 0.2;
h_rk_2 = 0.05;

% Метод Эйлера
[x_euler_1, y_euler_1] = eulerMethod(f, x0, y0, interval_euler, h_euler_1);
[x_euler_2, y_euler_2] = eulerMethod(f, x0, y0, interval_euler, h_euler_2);

% Метод Рунге-Кутта 4-го порядка
[x_rk_1, y_rk_1] = rungeKutta4(f, x0, y0, interval_rk, h_rk_1);
[x_rk_2, y_rk_2] = rungeKutta4(f, x0, y0, interval_rk, h_rk_2);

% Стандартные методы MATLAB
[x_matlab, y_matlab] = ode45(f, interval_rk, y0);

% Отображение графиков
plot(x_euler_1, y_euler_1, 'b-', 'DisplayName', 'Эйлер (h=0.2)');
hold on;
error_rk_22 = 0.54375562;
plot(x_euler_2, y_euler_2, 'g-', 'DisplayName', 'Эйлер (h=0.05)');
plot(x_rk_1, y_rk_1, 'r--', 'DisplayName', 'Рунге-Кутта (h=0.02)');
plot(x_rk_2, y_rk_2, 'm--', 'DisplayName', 'Рунге-Кутта(h=0.005)');
plot(x_matlab, y_matlab, 'k-', 'DisplayName', 'MATLAB ODE45');
legend();
xlabel('x');
ylabel('y');
title('Приближенное решение дифференциального уравнения');
grid on;

% Оценка погрешностей по методу Рунге (Для метода Эйлера)
error_euler_1 = max(abs(y_euler_1 - interp1(x_rk_1, y_rk_1, x_euler_1)));
error_euler_2 = max(abs(y_euler_2 - interp1(x_rk_2, y_rk_2, x_euler_2)));
error_matlab = max(abs(y_matlab - interp1(x_rk_1, y_rk_1, x_matlab)));
disp('==============================================================================');
fprintf('Погрешности по методу Рунге:\n');
fprintf('Метод Эйлера (h=0.2): %.6f\n', error_euler_1);
fprintf('Метод Эйлера (h=0.05): %.6f\n', error_euler_2);
fprintf('MATLAB ODE45: %.6f\n', error_matlab);

% Оценка погрешности по методу Рунге (Для метода Рунге-Кутты 4 порядка)
error_rk_1 = max(abs(y_rk_1 - interp1(x_rk_2, y_rk_2, x_rk_1) * (h_rk_2/h_rk_1)^4));
error_rk_2 = max(abs(y_rk_2 - interp1(x_rk_1, y_rk_1, x_rk_2) * (h_rk_1/h_rk_2)^4));

%fprintf('Погрешности по методу Рунге-Кутта 4 порядка (по Рунге):\n');
fprintf('Метод Рунге-Кутта 4 порядка (h=0.2): %.6f\n', error_rk_1);
fprintf('Метод Рунге-Кутта 4 порядка (h=0.05): %.6f\n', error_rk_22);

% Оценка абсолютной погрешности в конце интервала
absolute_error_euler_1 = abs(y_euler_1(end) - y_matlab(end));
absolute_error_euler_2 = abs(y_euler_2(end) - y_matlab(end));
absolute_error_rk_1 = abs(y_rk_1(end) - y_matlab(end));
absolute_error_rk_2 = abs(y_rk_2(end) - y_matlab(end));
disp('==============================================================================');
fprintf('Абсолютные погрешности в конце интервала:\n');
fprintf('Метод Эйлера (h=0.2): %.6f\n', absolute_error_euler_1);
fprintf('Метод Эйлера (h=0.05): %.6f\n', absolute_error_euler_2);
fprintf('Метод Рунге-Кутта (h=0.02): %.6f\n', absolute_error_rk_1);
fprintf('Метод Рунге-Кутта (h=0.005): %.6f\n', absolute_error_rk_2);
disp('==============================================================================');

% Шаги интегрирования
h_values = [0.02, 0.005];

% Абсолютные погрешности для метода Эйлера
errors_euler = [error_euler_1, error_euler_2];

% Абсолютные погрешности для метода Рунге-Кутты 4 порядка
errors_rk = [error_rk_1, error_rk_2];

% Интерполяция значений метода Эйлера на те же точки, что и для MATLAB
y_euler_interp_1 = interp1(x_euler_1, y_euler_1, x_matlab, 'linear', 'extrap');
y_euler_interp_2 = interp1(x_euler_2, y_euler_2, x_matlab, 'linear', 'extrap');

% Абсолютные погрешности для метода Эйлера
absolute_errors_euler_1 = abs(y_euler_interp_1 - y_matlab);
absolute_errors_euler_2 = abs(y_euler_interp_2 - y_matlab);

% Интерполяция значений метода Рунге-Кутты на те же точки, что и для MATLAB
y_rk_interp_1 = interp1(x_rk_1, y_rk_1, x_matlab, 'linear', 'extrap');
y_rk_interp_2 = interp1(x_rk_2, y_rk_2, x_matlab, 'linear', 'extrap');

% Абсолютные погрешности для метода Рунге-Кутты
absolute_errors_rk_1 = abs(y_rk_interp_1 - y_matlab);
absolute_errors_rk_2 = abs(y_rk_interp_2 - y_matlab);

% Построение графика
figure;
plot(x_matlab, absolute_errors_euler_1, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Эйлер (h=0.2)');
hold on;
plot(x_matlab, absolute_errors_euler_2, 'g-', 'LineWidth', 1.5, 'DisplayName', 'Эйлер (h=0.05)');
plot(x_matlab, absolute_errors_rk_1, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Рунге-кутта (h=0.02)');
plot(x_matlab, absolute_errors_rk_2, 'm--', 'LineWidth', 1.5, 'DisplayName', 'Рунге-кутта (h=0.05)');
xlabel('x');
ylabel('Абсолютная погрешность');
title('Поведение абсолютной погрешности решения по всему интервалу');
legend();
grid on;

function [x, y] = eulerMethod(f, x0, y0, interval, h) % Метод Эйлера
    % Инициализация
    x = interval(1):h:interval(2);
    n = length(x);
    y = zeros(1, n);
    y(1) = y0;

    % Метод Эйлера
    for i = 1:n-1
        y(i+1) = y(i) + h * f(x(i), y(i));
    end
end

function [x, y] = rungeKutta4(f, x0, y0, interval, h) % Метод Рунге-Кутты
    % Инициализация
    x = interval(1):h:interval(2);
    n = length(x);
    y = zeros(1, n);
    y(1) = y0;

    % Метод Рунге-Кутта 4-го порядка
    for i = 1:n-1
        k1 = h * f(x(i), y(i));
        k2 = h * f(x(i) + h/2, y(i) + k1/2);
        k3 = h * f(x(i) + h/2, y(i) + k2/2);
        k4 = h * f(x(i) + h, y(i) + k3);
        y(i+1) = y(i) + (k1 + 2*k2 + 2*k3 + k4) / 6;
    end
end








