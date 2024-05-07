clc; clear;
% Задание 2
f = @(x) sin(x(1) + x(2)) + (x(1) - x(2)).^2 - 1.5*x(1) + 2.5*x(2) + 1; % Функция для минимизации

% Градиент функции f
grad_f = @(x) [cos(x(1) + x(2)) + 2*(x(1) - x(2)) - 1.5; 
               cos(x(1) + x(2)) - 2*(x(1) - x(2)) + 2.5];

x0 = [0; 0]; % Начальное приближение

% Параметры метода
alpha = 0.01;
tol = 1e-4;
max_iter = 1000;

% Метод градиентного спуска
[xmin_custom, fmin_custom, steps_custom] = gradient_descent(f, grad_f, x0, alpha, tol, max_iter);

% Стандартный метод MATLAB
xmin_builtin = fminsearch(f, x0);
fmin_builtin = f(xmin_builtin);

% Вывод результатов
disp(' ')
disp('----------------------------------------------------------------------------------')
disp('Задание 2')
disp('Метод градиентного спуска:');
disp(['Минимум: ', num2str(fmin_custom)]);
disp(['Точка минимума: (', num2str(xmin_custom(1)), ', ', num2str(xmin_custom(2)), ')']);
disp('Стандартный метод MATLAB:');
disp(['Минимум: ', num2str(fmin_builtin)]);
disp(['Точка минимума: (', num2str(xmin_builtin(1)), ', ', num2str(xmin_builtin(2)), ')']);
disp('----------------------------------------------------------------------------------')

% График значений функции на отдельных шагах поиска экстремума
figure;
plot(1:length(steps_custom), steps_custom, 'bo-', 'LineWidth', 1.5);
hold on;
xlabel('Шаги');
ylabel('f(x,y)');
title('Значения функции на шагах градиентного спуска');
grid on;

% График сам найденного экстремума (собственный метод)
plot(length(steps_custom), fmin_custom, 'r*', 'LineWidth', 2);

% График сам найденного экстремума (fminsearch)
plot(length(steps_custom), f(xmin_builtin), 'g*', 'LineWidth', 2);

legend('Значения функции на шагах', 'Найденный минимум (метод градиентного спуска)', 'Найденный минимум (Стандартный метод MATLAB)');
hold off;

function [xmin, fmin, steps] = gradient_descent(f, grad_f, x0, alpha, tol, max_iter) % Метод градиентного спуска
    % Градиентный спуск для минимизации функции f
    % f - функция для минимизации
    % grad_f - градиент функции f
    % x0 - начальное приближение
    % alpha - скорость обучения
    % tol - допустимая ошибка
    % max_iter - максимальное количество итераций
    
    x = x0;
    steps = [];
    
    for iter = 1:max_iter
        grad = grad_f(x);
        x = x - alpha * grad;
        steps = [steps; f(x)];
        
        % Проверка на сходимость
        if norm(grad) < tol
            break;
        end
    end
    
    xmin = x;
    fmin = f(x);
end