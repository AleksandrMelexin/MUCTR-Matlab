clc; clear;
X = [4865 5065 5265 5465 5665 5865 6065 6265 6465 6665 6865 7065 7265 7465 7665 7865 8065 8265 8465 8665 8865];
Y = [2.799505654 2.865976811 2.936424531 3.010995579 3.089799019 3.172890855 3.260255665 3.351785366 3.44725572 3.546301795 3.648394432 3.752820666 3.858671918 3.964844364 4.070055767 4.172881981 4.27181405 4.365333422 4.451998938 4.530535861 4.599915455];

x_interp = 5500; % точка, в которой мы хотим вычислить производные

% Интерполяция
y_interp = interp1(X, Y, x_interp, 'spline');

% Построение графика интерполированной функции
plot(X, Y, 'o', x_interp, y_interp, 'rx');
xlabel('X');
ylabel('Y');
title('Интерполированная функция');
legend('Узлы данных', 'Интерполированная точка', 'Location', 'best');
grid on;

% Вычисление первой производной (численно)
dy_dx = gradient(Y) ./ gradient(X);
dy_dx_interp = interp1(X, dy_dx, x_interp, 'spline');

% Вычисление второй производной (численно)
d2y_dx2 = gradient(dy_dx) ./ gradient(X);
d2y_dx2_interp = interp1(X, d2y_dx2, x_interp, 'spline');

% Вывод результатов
disp(['Значение в точке x_interp: ', num2str(y_interp)]);
disp(['Ошибка в точке x_interp для 1-ой производной: ', num2str(dy_dx_interp)]);
disp(['Ошибка в точке x_interp для 2-ой производной: ', num2str(d2y_dx2_interp)]);


% Определение функции
f = @(x) x.^2 + exp(x + 3);

% Задание диапазона для x
x = linspace(-5, 5, 1000);

% Вычисление точных значения первой и второй производных
f_prime_exact = 2*x + exp(x + 3);
f_double_prime_exact = 2 + exp(x + 3);

% Простая формула для первой производной
h = 0.001;
f_prime_simple = (f(x + h) - f(x)) / h;

% Многоточечная формула для первой производной
f_prime_multistep = (-f(x + 2*h) + 8*f(x + h) - 8*f(x - h) + f(x - 2*h)) / (12*h);

% Простая формула для второй производной
f_double_prime_simple = (f(x + h) - 2*f(x) + f(x - h)) / h^2;

% Многоточечная формула для второй производной
f_double_prime_multistep = (-f(x + 2*h) + 16*f(x + h) - 30*f(x) + 16*f(x - h) - f(x - 2*h)) / (12*h^2);

% Вычисление погрешностей
error_prime_simple = abs(f_prime_exact - f_prime_simple);
error_prime_multistep = abs(f_prime_exact - f_prime_multistep);
error_double_prime_simple = abs(f_double_prime_exact - f_double_prime_simple);
error_double_prime_multistep = abs(f_double_prime_exact - f_double_prime_multistep);

% Построение графиков погрешностей
figure;

subplot(2, 1, 1);
plot(x, error_prime_simple, x, error_prime_multistep);
title('Погрешности первой производной');
legend('Простая формула', 'Многоточечная формула');
xlabel('x');
ylabel('Погрешность');

subplot(2, 1, 2);
plot(x, error_double_prime_simple, x, error_double_prime_multistep);
title('Погрешности второй производной');
legend('Простая формула', 'Многоточечная формула');
xlabel('x');
ylabel('Погрешность');



% Задание диапазона для x
x = linspace(-5, 5, 1000);

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