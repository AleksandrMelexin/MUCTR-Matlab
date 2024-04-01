clear; clc;
% Заданные уравнения
eq1 = @(x, y) (2 * x.^2 + 2 * y.^2 - 1).^5 - 22 * x.^2 * y.^3;
eq2 = @(x) -0.05 * x.^2 - 2.5 * x + 0.25;

% Начальное приближение для графического поиска
x_vals = linspace(-2, 2, 400);
y_vals = eq2(x_vals);

% Построение графиков
figure;
plot(x_vals, y_vals);
hold on;
fimplicit(eq1, [-2, 2, -2, 2], 'b');
xlabel('x');
ylabel('y');
title('Графики функций');
legend('y = -0.05 * x^2 - 2.5 * x + 0.25', '(2 * x^2 + 2 * y^2 - 1)^5 - 22 * x^2 * y^3 = 0');
grid on;

% Начальная точка для итерационных методов (взята из графика)
x0 = -1;
y0 = eq2(x0);

% Метод простых итераций
max_iter = 1000;
tolerance = 1e-6;

% Настройка параметров для метода Ньютона
options = optimset('Display', 'off');

% Решение системы уравнений методом простых итераций
[x_pi, y_pi, iter_pi] = simple_iteration_method(eq1, x0, y0, tolerance, max_iter);
if isnan(x_pi) || isnan(y_pi)
    disp('Метод простых итераций не сошелся');
else
    disp('Метод простых итераций:');
    disp(['x = ', num2str(x_pi)]);
    disp(['y = ', num2str(y_pi)]);
    disp(['Количество итераций: ', num2str(iter_pi)]);
end

% Решение системы уравнений методом Ньютона
[x_newton, y_newton, iter_newton] = newton_method(eq1, x0, y0, tolerance, max_iter, options);
if isnan(x_newton) || isnan(y_newton)
    disp('Метод Ньютона не сошелся');
else
    disp('Метод Ньютона:');
    disp(['x = ', num2str(x_newton)]);
    disp(['y = ', num2str(y_newton)]);
    disp(['Количество итераций: ', num2str(iter_newton)]);
end

% Решение с использованием встроенной функции fsolve
eq_system = @(vars) [eq1(vars(1), vars(2)); eq2(vars(1)) - vars(2)];
[x_fsolve, ~, exitflag] = fsolve(eq_system, [x0, y0], options);
if exitflag <= 0
    disp('Функция fsolve не смогла найти решение');
else
    disp('Функция fsolve:');
    disp(['x = ', num2str(x_fsolve(1))]);
    disp(['y = ', num2str(x_fsolve(2))]);
end

% Решение методом символьной математики vpasolve
syms x_sym y_sym;
eq1_sym = (2 * x_sym^2 + 2 * y_sym^2 - 1)^5 - 22 * x_sym^2 * y_sym^3 == 0;
eq2_sym = y_sym == -0.05 * x_sym^2 - 2.5 * x_sym + 0.25;
[solX, solY] = vpasolve([eq1_sym, eq2_sym], [x_sym, y_sym]);
disp('Метод символьной математики vpasolve:');
disp('x = ');
disp(solX);
disp('y = ');
disp(solY);

% Построение графика с найденными решениями
scatter(double(solX), double(solY), 'c', 'filled');
legend('y = -0.05 * x^2 - 2.5 * x + 0.25', '(2 * x^2 + 2 * y^2 - 1)^5 - 22 * x^2 * y^3 = 0', ...
    'vpasolve');
hold off;

% Определение функции для метода простых итераций
function [x, y, iter] = simple_iteration_method(eq1, x0, y0, tolerance, max_iter)
    iter = 0;
    while iter < max_iter
        x = eq1(x0, y0);
        y = -0.05 * x0.^2 - 2.5 * x0 + 0.25;
        if norm([x, y] - [x0, y0]) < tolerance
            return;
        end
        x0 = x;
        y0 = y;
        iter = iter + 1;
    end
    x = NaN;
    y = NaN;
end

% Определение функции для метода Ньютона
function [x, y, iter] = newton_method(eq1, x0, y0, tolerance, max_iter, options)
    iter = 0;
    while iter < max_iter
        J = [10 * (2 * x0^2 + 2 * y0^2 - 1)^4 * 4 * x0 * (2 * x0^2 + 2 * y0^2 - 1), ...
            10 * (2 * x0^2 + 2 * y0^2 - 1)^4 * 4 * y0^2 - 66 * x0 * y0^3;
            -0.1 * x0 - 2.5, -1];
        F = [-22 * x0^2 * y0^3 + (2 * x0^2 + 2 * y0^2 - 1)^5;
            -0.05 * x0^2 - 2.5 * x0 + 0.25 - y0];
        delta = J \ (-F);
        x = x0 + delta(1);
        y = y0 + delta(2);
        if norm([delta(1), delta(2)]) < tolerance
            return;
        end
        x0 = x;
        y0 = y;
        iter = iter + 1;
    end
    x = NaN;
    y = NaN;
end