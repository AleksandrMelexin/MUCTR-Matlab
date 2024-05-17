clc; clear;
% Определяем функцию и ее аналитические производные
f = @(x) x.^2 + exp(x + 3); % Исходная функция
f_prime = @(x) 2*x + exp(x + 3); % Первая производная функции
f_double_prime = @(x) 2 + exp(x + 3); % Вторая производная функции

% Задаем шаг и точки для вычислений
h = 1e-3; % Увеличенный шаг для улучшения видимости ошибок
x = -5:0.1:35;

% Вычисление первой производной простой конечной разностью
f_prime_simple = @(x, h) (f(x + h) - f(x)) / h;

% Вычисление первой производной многоточечной конечной разностью (центральная разность)
f_prime_central = @(x, h) (f(x + h) - f(x - h)) / (2*h);

% Вычисление второй производной простой конечной разностью
f_double_prime_simple = @(x, h) (f(x + h) - 2*f(x) + f(x - h)) / h^2;

% Вычисление второй производной многоточечной конечной разностью (центральная разность)
f_double_prime_central = @(x, h) (f(x + h) - 2*f(x) + f(x - h)) / h^2;

% Вычисление ошибок
error_prime_simple = abs(f_prime(x) - arrayfun(@(x) f_prime_simple(x, h), x));
error_prime_central = abs(f_prime(x) - arrayfun(@(x) f_prime_central(x, h), x));
error_double_prime_simple = abs(f_double_prime(x) - arrayfun(@(x) f_double_prime_simple(x, h), x));
error_double_prime_central = abs(f_double_prime(x) - arrayfun(@(x) f_double_prime_central(x, h), x));

% Ограничение интервала для построения графиков
x_interval = 10:0.1:30;

% Индексы для ограничения графиков
idx = find(x >= 10 & x <= 30);

% Построение графиков ошибок
figure;
subplot(2,1,1);
semilogy(x(idx), error_prime_simple(idx), 'r', x(idx), error_prime_central(idx), 'b');
title('Ошибка аппроксимации первой производной');
legend('Простая формула', 'Многоточечная формула');
xlabel('x');
ylabel('Погрешность');

subplot(2,1,2);
semilogy(x(idx), error_double_prime_simple(idx), '-o', x(idx), error_double_prime_central(idx), '-x');
title('Ошибка аппроксимации второй производной');
legend('Простая формула', 'Многоточечная формула');
xlabel('x');
ylabel('Погрешность');

% Настройки для лучшей видимости графиков
xlim([10 30]); % Ограничение интервала по оси x

% Массив значений шага h
h_values = logspace(-6, 0, 100);

% Пустые массивы для хранения погрешностей
errors_simple_prime = zeros(size(h_values));
errors_multistep_prime = zeros(size(h_values));
errors_simple_double_prime = zeros(size(h_values));
errors_multistep_double_prime = zeros(size(h_values));

% Вычисление погрешностей для разных значений шага
for i = 1:length(h_values)
    h = h_values(i);
    % Простая формула для первой производной
    f_prime_simple = (f(x + h) - f(x)) / h;
    % Многоточечная формула для первой производной
    f_prime_multistep = (-f(x + 2*h) + 8*f(x + h) - 8*f(x - h) + f(x - 2*h)) / (12*h);
    % Простая формула для второй производной
    f_double_prime_simple = (f(x + h) - 2*f(x) + f(x - h)) / h^2;
    % Многоточечная формула для второй производной
    f_double_prime_multistep = (-f(x + 2*h) + 16*f(x + h) - 30*f(x) + 16*f(x - h) - f(x - 2*h)) / (12*h^2);
    % Вычисление погрешностей
    errors_simple_prime(i) = max(abs(2*x + exp(x + 3) - f_prime_simple));
    errors_multistep_prime(i) = max(abs(2*x + exp(x + 3) - f_prime_multistep));
    errors_simple_double_prime(i) = max(abs(2 + exp(x + 3) - f_double_prime_simple));
    errors_multistep_double_prime(i) = max(abs(2 + exp(x + 3) - f_double_prime_multistep));
end

% Построение графиков зависимости погрешностей от значения шага
figure;

subplot(2, 1, 1);
loglog(h_values, errors_simple_prime, h_values, errors_multistep_prime);
title('Зависимость погрешности первой производной от шага');
legend('Простая формула', 'Многоточечная формула');
xlabel('Шаг (h)');
ylabel('Погрешность');

subplot(2, 1, 2);
loglog(h_values, errors_simple_double_prime, h_values, errors_multistep_double_prime);
title('Зависимость погрешности второй производной от шага');
legend('Простая формула', 'Многоточечная формула');
xlabel('Шаг (h)');
ylabel('Погрешность');