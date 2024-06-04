clc; clear;
% Начальные приближения
start = [1.25, 2.5];

% Решение системы уравнений с использованием fsolve
options = optimoptions('fsolve', 'Display', 'iter');
solution = fsolve(@system_eqs, start, options);

% Вывод решения fsolve
disp('Решение системы уравнений с использованием fsolve:');
disp(['x = ', num2str(solution(1))]);
disp(['y = ', num2str(solution(2))]);

% Начальные приближения для метода Ньютона
x0 = 1.25;
y0 = 2.5;
tol = 1e-6;
max_iter = 100;

% Проверка сходимости для метода Ньютона
if newton_convergence_check(x0, y0)
    % Решение методом Ньютона
    [x_newton, y_newton] = newton_method(x0, y0, tol, max_iter);

    % Вывод решения методом Ньютона
    disp('Решение методом Ньютона:');
    disp(['x = ', num2str(x_newton)]);
    disp(['y = ', num2str(y_newton)]);
else
    disp('Метод Ньютона не сходится для данных начальных приближений.');
end

% Начальные приближения для метода простых итераций
x0 = 1.25;
y0 = 2.5;
tol = 1e-6;
max_iter = 100;

% Проверка сходимости для метода простых итераций
if simple_iterations_convergence_check(x0, y0)
    % Решение методом простых итераций
    [x_iter, y_iter] = simple_iterations(x0, y0, tol, max_iter);

    % Вывод решения методом простых итераций
    disp('Решение методом простых итераций:');
    disp(['x = ', num2str(x_iter)]);
    disp(['y = ', num2str(y_iter)]);
else
    disp('Метод простых итераций не сходится для данных начальных приближений.');
end

% Построение графика
fimplicit(@(x, y) sin(x + 1) - y - 1.2, [-2, 2, -2, 2]);
hold on;
fimplicit(@(x, y) 2*x + cos(y) - 2, [-2, 2, -2, 2]);

% Отметка найденного решения
plot(solution(1), solution(2), 'ro', 'MarkerSize', 10, 'DisplayName', 'fsolve');

if exist('x_newton', 'var') && exist('y_newton', 'var')
    plot(x_newton, y_newton, 'go', 'MarkerSize', 10, 'DisplayName', 'Newton');
end

if exist('x_iter', 'var') && exist('y_iter', 'var')
    plot(x_iter, y_iter, 'bo', 'MarkerSize', 10, 'DisplayName', 'Simple Iterations');
end

% Настройка графика
legend('show');
xlabel('x');
ylabel('y');
title('График уравнений и найденные решения');
grid on;
hold off;

% Определение системы уравнений
function F = system_eqs(vars)
    x = vars(1);
    y = vars(2);
    F(1) = sin(x + 1) - y - 1.2;
    F(2) = 2*x + cos(y) - 2;
end

function [x, y, converged] = newton_method(x0, y0, tol, max_iter)
    converged = false;
    for iter = 1:max_iter
        % Вычисление значений функций и их производных
        F1 = sin(x0 + 1) - y0 - 1.2;
        F2 = 2*x0 + cos(y0) - 2;
        J11 = cos(x0 + 1);
        J12 = -1;
        J21 = 2;
        J22 = -sin(y0);
        
        % Якобиан
        J = [J11, J12; J21, J22];
        
        % Вектор функций
        F = [F1; F2];
        
        % Решение системы линейных уравнений
        delta = J \ (-F);
        
        % Обновление значений
        x0 = x0 + delta(1);
        y0 = y0 + delta(2);
        
        % Проверка условия остановки
        if norm(delta) < tol
            converged = true;
            break;
        end
    end
    x = x0;
    y = y0;
end

% Метод простых итераций
function [x, y, converged] = simple_iterations(x0, y0, tol, max_iter)
    converged = false;
    for iter = 1:max_iter
        % Обновление значений
        x_new = (2 - cos(y0)) / 2;
        y_new = sin(x0 + 1) - 1.2;
        
        % Проверка условия остановки
        if abs(x_new - x0) < tol && abs(y_new - y0) < tol
            converged = true;
            break;
        end
        
        % Обновление переменных
        x0 = x_new;
        y0 = y_new;
    end
    x = x0;
    y = y0;
end

% Проверка сходимости для метода Ньютона
function converged = newton_convergence_check(x0, y0)
    % Вычисление значений производных
    J11 = cos(x0 + 1);
    J12 = -1;
    J21 = 2;
    J22 = -sin(y0);

    % Определение Якобиана
    J = [J11, J12; J21, J22];

    % Проверка невырожденности Якобиана
    if det(J) ~= 0
        converged = true;
    else
        converged = false;
    end
end

% Проверка сходимости для метода простых итераций
function converged = simple_iterations_convergence_check(x0, y0)
    % Вычисление значений производных
    g1_x = 0; % Производная функции g1 по x
    g1_y = sin(y0); % Производная функции g1 по y
    g2_x = cos(x0 + 1); % Производная функции g2 по x
    g2_y = 0; % Производная функции g2 по y

    % Определение матрицы Якоби
    J = [g1_x, g1_y; g2_x, g2_y];

    % Проверка нормы матрицы Якоби
    if norm(J, 'fro') < 1
        converged = true;
    else
        converged = false;
    end
end