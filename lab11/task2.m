clc; clear;
% Задание 2
f = @(x) sin(x(1) + x(2)) + (x(1) - x(2)).^2 - 1.5*x(1) + 2.5*x(2) + 1; % Функция для минимизации

% Градиент функции f
grad_f = @(x) [cos(x(1) + x(2)) + 2*(x(1) - x(2)) - 1.5; 
               cos(x(1) + x(2)) - 2*(x(1) - x(2)) + 2.5];

x0 = [0; 1]; % Начальное приближение

% Параметры метода
alpha = 0.01;
tol = 1e-4;
max_iter = 1000;

% Метод градиентного спуска
[xmin_custom, fmin_custom, steps_f, steps1_x, steps2_x] = gradient_descent(f, grad_f, x0, alpha, tol, max_iter);

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

% Генерация сетки точек для построения графика
[X, Y] = meshgrid(-5:0.1:5, -5:0.1:5);
Z = sin(X + Y) + (X - Y).^2 - 1.5*X + 2.5*Y + 1;

% Построение 3D графика функции f
figure;
surf(X, Y, Z);
hold on;

% Настройка параметров графика
xlabel('X');
ylabel('Y');
zlabel('f(X, Y)');
title('График функции f(X, Y)');

% Построение точек оптимума
scatter3(steps1_x(:,1), steps2_x(:,1), steps_f(:,1), 'b', 'filled');
hold on;

% Построение точки экстремума
scatter3(xmin_custom(1), xmin_custom(2), fmin_custom, 'r', 'filled');

% Включение легенды
legend('f(X, Y)', 'Текущие точки оптимума', 'Точка экстремума', 'Location', 'best');


function [xmin, fmin, f_steps, x1_steps, x2_steps] = gradient_descent(f, grad_f, x0, alpha, tol, max_iter) % Метод градиентного спуска
    % Градиентный спуск для минимизации функции f
    % f - функция для минимизации
    % grad_f - градиент функции f
    % x0 - начальное приближение
    % alpha - скорость обучения
    % tol - допустимая ошибка
    % max_iter - максимальное количество итераций
    
    x = x0;
    f_steps = [];
    x1_steps = [];
    x2_steps = [];
    
    for iter = 1:max_iter
        grad = grad_f(x);
        x = x - alpha * grad;
        f_steps = [f_steps; f(x)];
        x1_steps = [x1_steps; x(1)];
        x2_steps = [x2_steps; x(2)];
        
        % Проверка на сходимость
        if norm(grad) < tol
            break;
        end
    end
    xmin = x;
    fmin = f(x);
end