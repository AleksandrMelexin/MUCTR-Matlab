% Даны точки данных
mass_heatCapacity = [2.799505654 2.865976811 2.936424531 3.010995579 3.089799019 3.172890855 3.260255665 3.351785366 3.44725572 3.546301795 3.648394432 3.752820666 3.858671918 3.964844364 4.070055767 4.172881981 4.27181405 4.365333422 4.451998938 4.530535861 4.599915455];
pressure = [4865 5065 5265 5465 5665 5865 6065 6265 6465 6665 6865 7065 7265 7465 7665 7865 8065 8265 8465 8665 8865];
weights = [0.2 0.7 0.5 0.3 1.0 0.2 0.0 0.7 0.8 0.6 0.0 0.1 0.0 0.2 0.3 0.6];

% Определение степени полинома
n = numel(mass_heatCapacity);
degree = n - 2;

% Конечные разности
delta_mass_heatCapacity = diff(mass_heatCapacity);
delta_pressure = diff(pressure);
delta_weights = diff(weights);

% Объединение конечных разностей в один массив
all_deltas = [delta_mass_heatCapacity(:); delta_pressure(:); delta_weights(:)];

% Определение степени полинома
degree_poly = find(min(arrayfun(@(x) nnz(diff(sign(x))), all_deltas)) == ...
                    nnz(diff(sign(delta_mass_heatCapacity))));

disp(['Степень полинома: ', num2str(degree_poly)]);

% Построение аппроксимационного полинома без учета весовых коэффициентов с использованием определителя Вандермонда
A = vander(pressure);
coefficients_vander = A(:, end-degree_poly:end) \ mass_heatCapacity.';
polynomial_vander = polyval(coefficients_vander, pressure);

% Построение аппроксимационного полинома без учета весовых коэффициентов с использованием стандартных операторов MATLAB (polyfit)
coefficients_polyfit = polyfit(pressure, mass_heatCapacity, degree_poly);
polynomial_polyfit = polyval(coefficients_polyfit, pressure);

% Построение аппроксимационного полинома с учетом весовых коэффициентов с использованием функции spap2
polynomial_spap2 = spap2(degree_poly, 3, pressure, mass_heatCapacity, weights);

% Оценка точности аппроксимации
error_vander = norm(mass_heatCapacity - polynomial_vander, 2);
error_polyfit = norm(mass_heatCapacity - polynomial_polyfit, 2);
error_spap2 = norm(mass_heatCapacity - fnval(polynomial_spap2, pressure), 2);

disp(['Среднеквадратичная ошибка (определитель Вандермонда): ', num2str(error_vander)]);
disp(['Среднеквадратичная ошибка (polyfit): ', num2str(error_polyfit)]);
disp(['Среднеквадратичная ошибка (spap2): ', num2str(error_spap2)]);

% Построение графика
plot(pressure, mass_heatCapacity, 'b*', 'DisplayName', 'Исходные данные');
hold on;
plot(pressure, polynomial_vander, 'r', 'DisplayName', 'Аппроксимация (определитель Вандермонда)');
plot(pressure, polynomial_polyfit, 'g', 'DisplayName', 'Аппроксимация (polyfit)');
fnplt(polynomial_spap2, 'm', 'DisplayName', 'Аппроксимация (spap2)');
legend;
xlabel('Давление');
ylabel('Теплоёмкость массы');
title('Аппроксимация данных');
grid on;
hold off;