clc; clear;
% Данные
X = [4865 5065 5265 5465 5665 5865 6065 6265 6465 6665 6865 7065 7265 7465 7665 7865 8065 8265 8465 8665 8865];
Y = [2.799505654 2.865976811 2.936424531 3.010995579 3.089799019 3.172890855 3.260255665 3.351785366 3.44725572 3.546301795 3.648394432 3.752820666 3.858671918 3.964844364 4.070055767 4.172881981 4.27181405 4.365333422 4.451998938 4.530535861 4.599915455];

% Определяем точки для интерполяции
points_to_interpolate = [4965, 7565];

% Интерполяция сплайном
pp = spline(X, Y);

% Вычисление производных в точках
first_derivative = ppval(fnder(pp, 1), points_to_interpolate);
second_derivative = ppval(fnder(pp, 2), points_to_interpolate);

% Погрешность по формуле для первой производной
n = length(X) - 1;
h = X(2) - X(1);
delta_2 = diff(Y, 2);
delta_3 = diff(Y, 3);
f_3 = delta_3 / h^3;

% Погрешность первой производной в точке 4965 (между узлами 1 и 2)
i_1 = 1;
error_first_derivative_4965 = ((-1)^(n-i_1)) * ((factorial(i_1) * factorial(n-1)) / factorial(n+1)) * (h^n) * f_3(i_1);

% Погрешность первой производной в точке 7565 (между узлами 14 и 15)
i_2 = 14;
error_first_derivative_7565 = ((-1)^(n-i_2)) * ((factorial(i_2) * factorial(n-1)) / factorial(n+1)) * (h^n) * f_3(i_2-1);

% Погрешность по формуле для второй производной
n = 2; % порядок конечной разности 2
f_3 = delta_3 / h^3;

% Погрешность второй производной в точке 4965 (между узлами 1 и 2)
error_second_derivative_4965 = ((-1)^(n-1)) * ((factorial(1) * factorial(n-1)) / factorial(n+1)) * (h^n) * f_3(i_1);

% Погрешность второй производной в точке 7565 (между узлами 14 и 15)
error_second_derivative_7565 = ((-1)^(n-1)) * ((factorial(1) * factorial(n-1)) / factorial(n+1)) * (h^n) * f_3(i_2-1);

% Вывод результатов
disp('Первая производная:');
disp(first_derivative);
disp('Погрешности первой производной:');
disp([abs(error_first_derivative_4965), abs(error_first_derivative_7565)]);

disp('Вторая производная:');
disp(second_derivative);
disp('Погрешности второй производной:');
disp([abs(error_second_derivative_4965), abs(error_second_derivative_7565)]);

% Построение графиков
xx = linspace(min(X), max(X), 1000);
yy = ppval(pp, xx);
yy_first_derivative = ppval(fnder(pp, 1), xx);
yy_second_derivative = ppval(fnder(pp, 2), xx);

figure;
subplot(3,1,1);
plot(X, Y, 'o', xx, yy, '-');
title('Функция Y(X)');
xlabel('X');
ylabel('Y');
legend('Данные', 'Интерполяция', 'Location', 'Best');

subplot(3,1,2);
plot(xx, yy_first_derivative, '-');
title('Первая производная');
xlabel('X');
ylabel('dY/dX');

subplot(3,1,3);
plot(xx, yy_second_derivative, '-');
title('Вторая производная');
xlabel('X');
ylabel('d^2Y/dX^2');
