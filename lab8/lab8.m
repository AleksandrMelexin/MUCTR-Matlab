clear; clc;
mass_heatCapacity = [2.799505654 2.865976811 2.936424531 3.010995579 3.089799019 3.172890855 3.260255665 3.351785366 3.44725572 3.546301795 3.648394432 3.752820666 3.858671918 3.964844364 4.070055767 4.172881981 4.27181405 4.365333422 4.451998938 4.530535861 4.599915455];
pressure = [4865 5065 5265 5465 5665 5865 6065 6265 6465 6665 6865 7065 7265 7465 7665 7865 8065 8265 8465 8665 8865];
weights = [0.5 0.2 0.9 0.5 0.8 0.8 0.2 0.5 0.4 0.2 0.6 0.5 0.2 0.7 0.1 0.8 0.6 0.9 0.6 0.9 0.4];
% 1 
% Построение таблицы конечных разностей
dt = diff(mass_heatCapacity);
disp(['Таблица конечных разностей: ',dt]);
disp(dt);
diff_table = diffTable(pressure, mass_heatCapacity, transpose(weights));  % Вызов функции diffTable для вычисления таблицы конечных разностей

% Максимальная степень полинома (n-2), где n - количество узловых точек
max_degree = length(pressure) - 2;

% Определение степени полинома на основе максимальной разности в таблице
degree = find(max(abs(diff_table)) > 1e-6, 1, 'last');
degree = min(degree, max_degree);  % Ограничиваем степень полинома максимальным значением

% Вывод степени полинома
disp(['Рекомендуемая степень полинома: ', num2str(degree)]);

% 2 
% Построение аппроксимационного полинома
% без учёта весовых коэффициентов с использованием 
% определителя Вандермонда

% Создание матрицы Вандермонда
W = vander(pressure);
W = W(:,  end-4:end); % Выбор последних столбцов, соответствующих степеням полинома

% Решение системы линейных уравнений
coefficients = W \ mass_heatCapacity';

% Вычисление значений полинома для построения графика
p_density_w = polyval(flip(coefficients'), pressure);

% Построение графика
figure;
plot(pressure, mass_heatCapacity, 'o', pressure, p_density_w, '-');
xlabel('Давление в кПа');
ylabel('Теплоёмкость');
title(['Аппроксимационный полином(Вандермонд) степени ', num2str(degree)]);
legend('Экспериментальные данные', 'Аппроксимационный полином (Вандермонд)');
grid on;

% 3 
% Построить аппроксимационный полином 
% без учёта весовых коэффициентов с использованием 
% стандартных операторов MATLAB

% Аппроксимация полиномом
p = polyfit(pressure, mass_heatCapacity, degree);
% Вычисление значений полинома для построения графика
p_density_s = polyval(p, pressure);
% Построение графика
figure;
plot(pressure, mass_heatCapacity, 'o', pressure, p_density_s, 'r-');
xlabel('Давление в кПа');
ylabel('Теплоёмкость');
title(['Аппроксимационный полином степени ', num2str(degree)]);
legend('Экспериментальные данные', 'Аппроксимационный полином стандарт');
grid on;

% 4 
% Построить аппроксимационный полином 
% с учётом весовых коэффициентов с использованием 
% функции spap2

% Создание объекта сплайна
sp = spap2(pressure, 1, pressure, mass_heatCapacity, weights);

% Оценка значений плотности для построения графика
p_density2 = fnval(sp, pressure);

% Построение графика
figure;
plot(pressure, mass_heatCapacity, 'o', pressure, p_density2, 'g-');
xlabel('Давление в кПа');
ylabel('Теплоёмкость');
title(['Аппроксимационный полином(spap2) степени ', num2str(degree)]);
legend('Экспериментальные данные', 'Аппроксимационный полином spap2');
grid on;

% 5 Построить аппроксимационный полином 
% с учётом весовых коэффициентов с использованием 
% функции fminsearch

% Функция для минимизации (сумма квадратов разностей)
fun = @(coefficients) sum(weights .* (polyval(fliplr(coefficients'), pressure) - mass_heatCapacity).^2);

% Начальное приближение для коэффициентов полинома
initial_guess = rand(1, degree+1); % Случайное начальное приближение

% Минимизация функции с использованием fminsearch
optimal_coefficients = fminsearch(fun, initial_guess);

% Вычисление значений полинома для построения графика
p_density_f = polyval(flip(optimal_coefficients'), pressure);

% Построение графика
figure;
plot(pressure, mass_heatCapacity, 'o', pressure, p_density_f, 'b-');
xlabel('Давление в кПа');
ylabel('Теплоёмкость');
title(['Аппроксимационный полином(fminsearch) степени ', num2str(degree)]);
legend('Экспериментальные данные', 'Аппроксимационный полином fminsearch');
grid on;

% 6 
% Оценить точность аппроксимации
% Вычисление MSE для каждого метода
mse_polyfit = mean((polyval(p, pressure) - mass_heatCapacity).^2);
mse_vander = mean((polyval(flip(coefficients'), pressure) - mass_heatCapacity).^2);
mse_spap2 = mean((fnval(sp, pressure) - mass_heatCapacity).^2);
mse_fminsearch = mean((polyval(flip(optimal_coefficients'), pressure) - mass_heatCapacity).^2);

% Вывод результатов
disp('Среднеквадратичная ошибка (MSE) для каждого метода:');
disp(['1. Полином методом наименьших квадратов (без учета весовых коэффициентов): ', num2str(mse_polyfit)]);
disp(['2. Полином с использованием определителя Вандермонда: ', num2str(mse_vander)]);
disp(['3. Сплайны с использованием функции spap2: ', num2str(mse_spap2)]);
disp(['4. Аппроксимация с использованием функции fminsearch: ', num2str(mse_fminsearch)]);

% 1 и 2, 14 и 15 из 3 пункта
% Вычисление абсолютной погрешности для методов 1 и 2
abs_error_1_2_polyfit = abs(polyval(p, pressure(2)) - mass_heatCapacity(2));
abs_error_1_2_vander = abs(polyval(flip(coefficients'), pressure(2)) - mass_heatCapacity(2));

% Вычисление абсолютной погрешности для методов 14 и 15
abs_error_14_15_polyfit = abs(polyval(p, pressure(15)) - mass_heatCapacity(15));
abs_error_14_15_vander = abs(polyval(flip(coefficients'), pressure(15)) - mass_heatCapacity(15));

% Вычисление среднего абсолютного значения погрешности для каждого метода
mean_abs_error_1_2_polyfit = mean(abs_error_1_2_polyfit);
mean_abs_error_1_2_vander = mean(abs_error_1_2_vander);
mean_abs_error_14_15_polyfit = mean(abs_error_14_15_polyfit);
mean_abs_error_14_15_vander = mean(abs_error_14_15_vander);

% Вывод результатов
disp('Среднее абсолютное значение погрешности для каждого метода:');
disp(['1. Полином методом наименьших квадратов (без учета весовых коэффициентов), между 1 и 2 точками: ', num2str(mean_abs_error_1_2_polyfit)]);
disp(['2. Полином с использованием определителя Вандермонда, между 1 и 2 точками: ', num2str(mean_abs_error_1_2_vander)]);
disp(['3. Полином методом наименьших квадратов (без учета весовых коэффициентов), между 14 и 15 точками: ', num2str(mean_abs_error_14_15_polyfit)]);
disp(['4. Полином с использованием определителя Вандермонда, между 14 и 15 точками: ', num2str(mean_abs_error_14_15_vander)]);

% 7 
% Построить на одном графике полученные функций 
% с вынесенными на них узловыми точками в виде звёздочек, разными цветами, 
% добавив легенду.

figure;
hold on
grid on
% точки и Вандермонд
plot(pressure, mass_heatCapacity, 'o', pressure, p_density_w, '-');
% стандарт
plot(pressure, p_density_s, '-');
% спап2
plot(pressure, p_density2, '-');
% фминсерч
plot(pressure, p_density_f, '-');
xlabel('Давление в кПа');
ylabel('Теплоёмкость');
title('Аппроксимационные полиномы');
legend('Экспериментальные данные', 'Полином + Вандермонд', 'стандартные методы', 'spap2', 'fminsearch');
hold off

% 8 
% Построить функцию, аппроксимирующую экспериментальные данные, 
% не в виде полинома, в виде другой функции 
% с учётом весовых коэффициентов с использованием функции fminsearch

% функция к минимизации(мнк)
obj_function = @(params) sum(weights .* (mass_heatCapacity - (params(1) .* pressure + params(2)) ./ (params(3) .* pressure + params(4))).^2);

% Догадка
initial_guess = [1, 1, 1, 1];

% оптимизация fminsearch
params = fminsearch(obj_function, initial_guess);

% Инициализация параметров
a = params(1);
b = params(2);
c = params(3);
d = params(4);

% Вывод подобранных параметров
disp(['Optimized Parameters: a = ', num2str(a), ', b = ', num2str(b), ', c = ', num2str(c), ', d = ', num2str(d)]);

% График данных и кривой
p_range = linspace(min(pressure), max(pressure), 100);
fitted_density = (a .* p_range + b) ./ (c .* p_range + d);

figure;
plot(pressure, mass_heatCapacity, 'bo', p_range, fitted_density, 'm-');
grid on
xlabel('Давление (кПа)');
ylabel('Теплоёмкость');
legend('Экспериментальные данный', 'Кривая');
title('Аппроксимационный полином без полинома');

% 9 
% Построить функцию, аппроксимирующую экспериментальные данные, 
% в виде полинома Чебышёва, с учётом весовых коэффициентов

% Определим степень полинома Чебышева
degree = 5;

% Вычисляем узлы 
x = 2 * ((pressure - min(pressure)) / (max(pressure) - min(pressure))) - 1;

% Вычисление многочлена с помощью нормализованных мнк
p_coefficients = polyfit(x, mass_heatCapacity, degree);

% многочлен в исходных точках данных
fitted_density = polyval(p_coefficients, x);

% Построение исходных данных и подогнанной кривой
figure;
plot(pressure, mass_heatCapacity, 'bo', pressure, fitted_density, 'c-');
grid on
xlabel('Давление (кПа)');
ylabel('Теплоёмкость');
legend('Экспериментальные данный', 'Полином');
title('Аппроксимационный полином Чебышева');


function diff_table = diffTable(x, y, w)
    % Инициализация таблицы конечных разностей
    n = length(x);
    diff_table = zeros(n, n);
    diff_table(:, 1) = w .* y';

    % Вычисление конечных разностей
    for j = 2:n
        for i = 1:n-j+1
            w_factor = w(i) / w(i+j-1);  % Коэффициент весов для текущей разности
            diff_table(i, j) = w_factor * (diff_table(i+1, j-1) - diff_table(i, j-1));
        end
    end
end
