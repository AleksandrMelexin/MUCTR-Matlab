clc; clear;
% Задание 1
% Начальный интервал
a = 10;
b = 30;
d0 = 5; % Начальное приближение
tol = 1e-4; % Точность и максимальное количество итераций
max_iter = 100;

% Вызов методов Ньютона, Золотого чесения и парабол
min_d_newton = newton_method(d0, tol, max_iter);
min_d_golden = golden_section_method(a, b, tol, max_iter);
min_d_parabolic = parabolic_method(a, b, tol, max_iter);

% Вызов функции MATLAB
options = optimset('Display', 'off'); 
[min_d_matlab, min_heat_loss] = fminsearch(@calculate_heat_loss, d0, options);

% Вывод результата
disp('----------------------------------------------------------------------------------')
disp('Задание 1')
disp(['Оптимальное значение d (метод Ньютона): ', num2str(min_d_newton)]);
disp(['Оптимальное значение d (метод золотого сечения): ', num2str(min_d_golden)]);
disp(['Оптимальное значение d (метод парабол): ', num2str(min_d_parabolic)]);
disp(['Оптимальное значение d (стандартная функция MATLAB): ', num2str(min_d_matlab)]);
disp('----------------------------------------------------------------------------------')

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

function heat_loss = calculate_heat_loss(d) % Функция для вычисления потерь тепла
    h = (4*1200)/(pi*d.^2); % выражение высоты исходя из общего объёма
    % Коэффициенты теплоотдачи
    h_bottom = 0.35; % Для дна
    h_wall = 1; % Для стенок
    h_lid = 0.68; % Для крышки
    
    % Площади поверхностей
    S_bottom = pi * (d^2) / 4;
    S_wall = pi * d * h;
    S_lid = S_bottom;
    
    % Тепловые потери через каждую поверхность
    Q_bottom = h_bottom * S_bottom;
    Q_wall = h_wall * S_wall;
    Q_lid = h_lid * S_lid;
    
    % Общие тепловые потери
    heat_loss = Q_bottom + Q_wall + Q_lid;
end

function min_d = newton_method(d0, tol, max_iter) % метод Ньютона
    d = d0;
    iter = 0;
    while iter < max_iter
        h_loss = calculate_heat_loss(d);
        h = 1e-4; % малое приращение для численного дифференцирования
        d_prime = (calculate_heat_loss(d + h) - calculate_heat_loss(d)) / h;
        d_double_prime = (calculate_heat_loss(d + h) - 2 * h_loss + calculate_heat_loss(d - h)) / h^2;
        d = d - d_prime / d_double_prime;
        if abs(d_prime) < tol
            break;
        end
        iter = iter + 1;
    end
    min_d = d;
end

function min_d = golden_section_method(a, b, tol, max_iter) % метод золотого сечения
    rho = (sqrt(5) - 1) / 2;
    d = a + rho * (b - a);
    c = b - rho * (b - a);
    iter = 0; % Инициализация переменной iter
    while abs(b - a) > tol && iter < max_iter
        if calculate_heat_loss(c) < calculate_heat_loss(d)
            b = d;
            d = c;
            c = b - rho * (b - a);
        else
            a = c;
            c = d;
            d = a + rho * (b - a);
        end
        iter = iter + 1; % Увеличение счетчика итераций
    end
    min_d = (a + b) / 2;
end

function min_d = parabolic_method(a, b, tol, max_iter) % метод парабол
    x = [a, (a + b) / 2, b];
    f = zeros(1, 3);
    f(1) = calculate_heat_loss(x(1));
    f(2) = calculate_heat_loss(x(2));
    f(3) = calculate_heat_loss(x(3));
    iter = 0;
    while abs(b - a) > tol && iter < max_iter
        A = ((x(2) - x(3)) * f(1) + (x(3) - x(1)) * f(2) + (x(1) - x(2)) * f(3)) / ...
            ((x(1) - x(2)) * (x(2) - x(3)) * (x(3) - x(1)));
        B = (f(2) - f(1)) / (x(2) - x(1)) - (x(1) + x(2)) * A;
        C = f(1) - A * x(1)^2 - B * x(1);
        min_x = -B / (2 * A);
        if min_x < x(2)
            if calculate_heat_loss(min_x) < f(2)
                x(3) = x(2);
                x(2) = min_x;
                f(3) = f(2);
                f(2) = calculate_heat_loss(min_x);
            else
                x(1) = min_x;
                f(1) = calculate_heat_loss(min_x);
            end
        else
            if calculate_heat_loss(min_x) < f(2)
                x(1) = x(2);
                x(2) = min_x;
                f(1) = f(2);
                f(2) = calculate_heat_loss(min_x);
            else
                x(3) = min_x;
                f(3) = calculate_heat_loss(min_x);
            end
        end
        iter = iter + 1;
    end
    min_d = x(2);
end

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