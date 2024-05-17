clc; clear;
% Начальные условия для первой системы
y0_system1 = [1; 1];
x_span_system1 = [0, 2];

% Начальные условия для второй системы
y0_system2 = [1; 1];
x_span_system2 = [0, 2];
    
% Шаг
h = 0.1;
    
% Вычисление решения первой системы стандартным методом MATLAB
ode_options = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);
[x_system1, y_system1] = ode45(@ode_func_system1, x_span_system1, y0_system1, ode_options);
y1_standard_system1 = y_system1(:, 1)';
y2_standard_system1 = y_system1(:, 2)';
    
% Решение первой системы методом Эйлера, Модифицированным методом Эйлера и методом Рунге-Кутты 4-го порядка
[x_euler_system1, y1_euler_system1, y2_euler_system1] = euler_method(y0_system1, x_span_system1, h);
[x_mod_euler_system1, y1_mod_euler_system1, y2_mod_euler_system1] = modified_euler_method(y0_system1, x_span_system1, h);
[x_rk4_system1, y1_rk4_system1, y2_rk4_system1] = rk4_method(y0_system1, x_span_system1, h);
    
% Вычисление числа жесткости для первой системы
stiffness_euler_system1 = compute_stiffness(x_euler_system1, y1_euler_system1, y2_euler_system1);
stiffness_mod_euler_system1 = compute_stiffness(x_mod_euler_system1, y1_mod_euler_system1, y2_mod_euler_system1);
stiffness_rk4_system1 = compute_stiffness(x_rk4_system1, y1_rk4_system1, y2_rk4_system1);
stiffness_standard_system1 = compute_stiffness(x_system1, y1_standard_system1, y2_standard_system1);
    
% Построение графика числа жесткости для первой системы
figure;
hold on;
plot(x_euler_system1, stiffness_euler_system1, 'b', 'LineWidth', 2);
plot(x_mod_euler_system1, stiffness_mod_euler_system1, 'r', 'LineWidth', 2);
plot(x_rk4_system1, stiffness_rk4_system1, 'g', 'LineWidth', 2);
plot(x_system1, stiffness_standard_system1, 'k', 'LineWidth', 2);
xlabel('x');
ylabel('Stiffness');
title('Число жесткости для первой системы');
legend('Эйлер', 'Модифицированный Эйлер', 'Рунге-Кутты', 'MATLAB');
grid on;
    
% Вычисление решения второй системы стандартным методом MATLAB
[x_system2, y_system2] = ode45(@ode_func_system2, x_span_system2, y0_system2, ode_options);
y1_standard_system2 = y_system2(:, 1)';
y2_standard_system2 = y_system2(:, 2)';
    
% Решение второй системы методом Эйлера, Модифицированным методом Эйлера и методом Рунге-Кутты 4-го порядка
[x_euler_system2, y1_euler_system2, y2_euler_system2] = euler_method(y0_system2, x_span_system2, h);
[x_mod_euler_system2, y1_mod_euler_system2, y2_mod_euler_system2] = modified_euler_method(y0_system2, x_span_system2, h);
[x_rk4_system2, y1_rk4_system2, y2_rk4_system2] = rk4_method(y0_system2, x_span_system2, h);
    
% Вычисление числа жесткости для второй системы
stiffness_euler_system2 = compute_stiffness(x_euler_system2, y1_euler_system2, y2_euler_system2);
stiffness_mod_euler_system2 = compute_stiffness(x_mod_euler_system2, y1_mod_euler_system2, y2_mod_euler_system2);
stiffness_rk4_system2 = compute_stiffness(x_rk4_system2, y1_rk4_system2, y2_rk4_system2);
stiffness_standard_system2 = compute_stiffness(x_system2, y1_standard_system2, y2_standard_system2);
    
% Построение графика числа жесткости для второй системы
figure;
hold on;
plot(x_euler_system2, stiffness_euler_system2, 'b', 'LineWidth', 2);
plot(x_mod_euler_system2, stiffness_mod_euler_system2, 'r', 'LineWidth', 2);
plot(x_rk4_system2, stiffness_rk4_system2, 'g', 'LineWidth', 2);
plot(x_system2, stiffness_standard_system2, 'k', 'LineWidth', 2);
xlabel('x');
ylabel('Stiffness');
title('Число жесткости для второй системы');
legend('Эйлер', 'Модифицированный Эйлер', 'Рунге-Кутты', 'MATLAB');
grid on;
    
% Построение графиков для первой системы
plot_results(x_system1, y1_standard_system1, y2_standard_system1, x_euler_system1, y1_euler_system1, y2_euler_system1, x_mod_euler_system1, y1_mod_euler_system1, y2_mod_euler_system1, x_rk4_system1, y1_rk4_system1, y2_rk4_system1, 'Решение первой системы');
    
% Построение графиков для второй системы
plot_results(x_system2, y1_standard_system2, y2_standard_system2, x_euler_system2, y1_euler_system2, y2_euler_system2, x_mod_euler_system2, y1_mod_euler_system2, y2_mod_euler_system2, x_rk4_system2, y1_rk4_system2, y2_rk4_system2, 'Решение второй системы');

function dydx = ode_func_system1(x, y)
    dydx = zeros(2, 1);
    dydx(1) = y(1)*exp(-x^2) + x*y(2);
    dydx(2) = 3*x - y(1) + 2*y(2);
end

function dydx = ode_func_system2(x, y)
    dydx = zeros(2, 1);
    dydx(1) = y(1)*exp(x^2) + x*y(2);
    dydx(2) = 3*x - y(1) + 2*y(2);
end

function plot_results(x_standard, y1_standard, y2_standard, x_euler, y1_euler, y2_euler, x_mod_euler, y1_mod_euler, y2_mod_euler, x_rk4, y1_rk4, y2_rk4, title_str)
    figure;
    hold on;
    plot(x_standard, y1_standard, 'k', 'LineWidth', 2);
    plot(x_standard, y2_standard, 'k--', 'LineWidth', 2);
    plot(x_euler, y1_euler, 'b', 'LineWidth', 2);
    plot(x_euler, y2_euler, 'b--', 'LineWidth', 2);
    plot(x_mod_euler, y1_mod_euler, 'r', 'LineWidth', 2);
    plot(x_mod_euler, y2_mod_euler, 'r--', 'LineWidth', 2);
    plot(x_rk4, y1_rk4, 'g', 'LineWidth', 2);
    plot(x_rk4, y2_rk4, 'g--', 'LineWidth', 2);
    xlabel('x');
    ylabel('y');
    title(title_str);
    legend('y1 - MATLAB', 'y2 - MATLAB', 'y1 - Эйлер', 'y2 - Эйлер', 'y1 - Модифицированный Эйлер', 'y2 - Модифицированный Эйлер', 'y1 - Рунге-Кутты', 'y2 - Рунге-Кутты');
    grid on;
end

function [x, y1, y2] = euler_method(y0, x_span, h)
    x = x_span(1):h:x_span(2);
    y1 = zeros(size(x));
    y2 = zeros(size(x));
    y1(1) = y0(1);
    y2(1) = y0(2);
    for i = 1:(length(x)-1)
        y1_prime = y1(i)*exp(-x(i)^2) + x(i)*y2(i);
        y2_prime = 3*x(i) - y1(i) + 2*y2(i);
        y1(i+1) = y1(i) + h * y1_prime;
        y2(i+1) = y2(i) + h * y2_prime;
    end
end

function [x, y1, y2] = modified_euler_method(y0, x_span, h)
    x = x_span(1):h:x_span(2);
    y1 = zeros(size(x));
    y2 = zeros(size(x));
    y1(1) = y0(1);
    y2(1) = y0(2);
    for i = 1:(length(x)-1)
        y1_prime = y1(i)*exp(-x(i)^2) + x(i)*y2(i);
        y2_prime = 3*x(i) - y1(i) + 2*y2(i);
        y1_temp = y1(i) + h * y1_prime;
        y2_temp = y2(i) + h * y2_prime;
        y1(i+1) = y1(i) + h/2 * (y1_prime + (y1_temp*exp(-(x(i)+h)^2) + (x(i)+h)*y2_temp));
        y2(i+1) = y2(i) + h/2 * (y2_prime + (3*(x(i)+h) - y1_temp + 2*y2_temp));
    end
end

function [x, y1, y2] = rk4_method(y0, x_span, h)
    x = x_span(1):h:x_span(2);
    y1 = zeros(size(x));
    y2 = zeros(size(x));
    y1(1) = y0(1);
    y2(1) = y0(2);
    for i = 1:(length(x)-1)
        k1_y1 = y1(i)*exp(-x(i)^2) + x(i)*y2(i);
        k1_y2 = 3*x(i) - y1(i) + 2*y2(i);
        
        k2_y1 = (y1(i) + h/2*k1_y1)*exp(-(x(i)+h/2)^2) + (x(i)+h/2)*(y2(i) + h/2*k1_y2);
        k2_y2 = 3*(x(i)+h/2) - (y1(i) + h/2*k1_y1) + 2*(y2(i) + h/2*k1_y2);
        
        k3_y1 = (y1(i) + h/2*k2_y1)*exp(-(x(i)+h/2)^2) + (x(i)+h/2)*(y2(i) + h/2*k2_y2);
        k3_y2 = 3*(x(i)+h/2) - (y1(i) + h/2*k2_y1) + 2*(y2(i) + h/2*k2_y2);
        
        k4_y1 = (y1(i) + h*k3_y1)*exp(-(x(i)+h)^2) + (x(i)+h)*(y2(i) + h*k3_y2);
        k4_y2 = 3*(x(i)+h) - (y1(i) + h*k3_y1) + 2*(y2(i) + h*k3_y2);
        
        y1(i+1) = y1(i) + h/6 * (k1_y1 + 2*k2_y1 + 2*k3_y1 + k4_y1);
        y2(i+1) = y2(i) + h/6 * (k1_y2 + 2*k2_y2 + 2*k3_y2 + k4_y2);
    end
end

function stiffness = compute_stiffness(x, y1, y2)
    stiffness = zeros(size(x));
    for i = 1:length(x)
        y1_prime = y1(i)*exp(-x(i)^2) + x(i)*y2(i);
        y2_prime = 3*x(i) - y1(i) + 2*y2(i);
        stiffness(i) = max(abs([y1_prime, y2_prime]));
    end
end
