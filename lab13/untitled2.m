clc; clear;
h = 0.1;
x_end_1 = 1;
x_end_2 = 3;
    
% Начальные условия
y0 = [1; 1];
    
% Первая система
fprintf('Решение первой системы:\n');
% Методы Эйлера, модифицированный Эйлера и Рунге-Кутта 4 порядка
[x_euler, y_euler] = euler(@first_system, y0, h, x_end_1);
[x_mod_euler, y_mod_euler] = modified_euler(@first_system, y0, h, x_end_1);
[x_rk4, y_rk4] = rk4(@first_system, y0, h, x_end_1);
% Стандартный метод MATLAB
[x_ode45, y_ode45] = ode45(@first_system, [0, x_end_1], y0);
    
% Построение графика
figure;
plot(x_euler, y_euler(1, :), 'b', x_mod_euler, y_mod_euler(1, :), 'r', x_rk4, y_rk4(1, :), 'g', x_ode45, y_ode45(:, 1), 'k');
legend('Эйлер', 'Модифицированный Эйлер', 'Рунге-Кутта 4 порядка', 'ODE45');
title('Первая система: y_1 от x');
xlabel('x');
ylabel('y_1');
    
figure;
plot(x_euler, y_euler(2, :), 'b', x_mod_euler, y_mod_euler(2, :), 'r', x_rk4, y_rk4(2, :), 'g', x_ode45, y_ode45(:, 2), 'k');
legend('Эйлер', 'Модифицированный Эйлер', 'Рунге-Кутта 4 порядка', 'ODE45');
title('Первая система: y_2 от x');
xlabel('x');
ylabel('y_2');
    
% Вторая система
fprintf('Решение второй системы:\n');
[x_euler_2, y_euler_2] = euler(@second_system, y0, h, x_end_2);
[x_implicit_euler, y_implicit_euler] = implicit_euler(@second_system, y0, h, x_end_2);
[x_ode45_2, y_ode45_2] = ode45(@second_system, [0, x_end_2], y0);
    
% Построение графика
figure;
plot(x_euler_2, y_euler_2(1, :), 'b', x_implicit_euler, y_implicit_euler(1, :), 'r', x_ode45_2, y_ode45_2(:, 1), 'k');
legend('Эйлер', 'Неявный Эйлер', 'ODE45');
title('Вторая система: y_1 от x');
xlabel('x');
ylabel('y_1');
    
figure;
plot(x_euler_2, y_euler_2(2, :), 'b', x_implicit_euler, y_implicit_euler(2, :), 'r', x_ode45_2, y_ode45_2(:, 2), 'k');
legend('Эйлер', 'Неявный Эйлер', 'ODE45');
title('Вторая система: y_2 от x');
xlabel('x');
ylabel('y_2');

function dy = first_system(x, y)
    dy = zeros(2, 1);
    dy(1) = y(1) * exp(-x^2) + x * y(2);
    dy(2) = 3 * x - y(1) + 2 * y(2);
end

function dy = second_system(x, y)
    dy = zeros(2, 1);
    dy(1) = y(1) * exp(x^2) + x * y(2);
    dy(2) = 3 * x - y(1) + 2 * y(2);
end

function [x, y] = euler(odefun, y0, h, x_end)
    n = floor(x_end / h);
    x = linspace(0, x_end, n + 1);
    y = zeros(length(y0), n + 1);
    y(:, 1) = y0;
    for i = 1:n
        y(:, i + 1) = y(:, i) + h * odefun(x(i), y(:, i));
    end
end

function [x, y] = modified_euler(odefun, y0, h, x_end)
    n = floor(x_end / h);
    x = linspace(0, x_end, n + 1);
    y = zeros(length(y0), n + 1);
    y(:, 1) = y0;
    for i = 1:n
        k1 = odefun(x(i), y(:, i));
        k2 = odefun(x(i) + h, y(:, i) + h * k1);
        y(:, i + 1) = y(:, i) + h / 2 * (k1 + k2);
    end
end

function [x, y] = rk4(odefun, y0, h, x_end)
    n = floor(x_end / h);
    x = linspace(0, x_end, n + 1);
    y = zeros(length(y0), n + 1);
    y(:, 1) = y0;
    for i = 1:n
        k1 = odefun(x(i), y(:, i));
        k2 = odefun(x(i) + h / 2, y(:, i) + h / 2 * k1);
        k3 = odefun(x(i) + h / 2, y(:, i) + h / 2 * k2);
        k4 = odefun(x(i) + h, y(:, i) + h * k3);
        y(:, i + 1) = y(:, i) + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4);
    end
end

function [x, y] = implicit_euler(odefun, y0, h, x_end)
    n = floor(x_end / h);
    x = linspace(0, x_end, n + 1);
    y = zeros(length(y0), n + 1);
    y(:, 1) = y0;
    for i = 1:n
        y(:, i + 1) = fsolve(@(Y) Y - y(:, i) - h * odefun(x(i + 1), Y), y(:, i));
    end
end