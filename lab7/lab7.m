clc; clear;
mass_heatCapacity = [2.799505654 2.865976811 2.936424531 3.010995579 3.089799019 3.172890855 3.260255665 3.351785366 3.44725572 3.546301795 3.648394432 3.752820666 3.858671918 3.964844364 4.070055767 4.172881981 4.27181405 4.365333422 4.451998938 4.530535861 4.599915455];
pressure = [4865 5065 5265 5465 5665 5865 6065 6265 6465 6665 6865 7065 7265 7465 7665 7865 8065 8265 8465 8665 8865];
weights = [0.2 0.7 0.5 0.3 1.0 0.2 0.0 0.7 0.8 0.6 0.0 0.1 0.0 0.2 0.3 0.6];

% Интерполяция данных
p = polyfit(pressure, mass_heatCapacity, length(pressure)-1); 
vandermonde_det = det(vander(pressure));
interp_temp = pressure;
interpolation = polyval(p, interp_temp);
[pl] = buildCanonicalPoly(mass_heatCapacity, pressure);

% Сплайн-интерполяция
s = spline(pressure, mass_heatCapacity, interp_temp); 

% Оценка погрешности интерполяции
midpoint1 = (pressure(14) + pressure(15)) / 2; 
midpoint2 = (pressure(1) + pressure(2)) / 2; 
interp_value1 = polyval(p, midpoint1); 
interp_value2 = polyval(p, midpoint2); 
actual_value1 = interp1(pressure, mass_heatCapacity, midpoint1); 
actual_value2 = interp1(pressure, mass_heatCapacity, midpoint2); 
error1 = abs(interp_value1 - actual_value1); 
error2 = abs(interp_value2 - actual_value2);

% построение таблицы конечных разностей
table_diffs = zeros(length(mass_heatCapacity));
table_diffs(:, 1) = mass_heatCapacity'; % Заполняем первый столбец таблицы значениями mass_heatCapacity
for j = 2:length(mass_heatCapacity)
    for i = 1:length(mass_heatCapacity)-j+1
        table_diffs(i, j) = table_diffs(i+1, j-1) - table_diffs(i, j-1); % Вычисляем конечные разности для каждого столбца
    end
end
disp('Таблица конечных разностей:');
disp(table_diffs);

% нахождение оптимальной степени полинома
fprintf('Массив разностей: ');
differences = zeros(1, size(table_diffs, 2));
for j = 1:size(table_diffs, 2)
    column = table_diffs(:, j);
    differences(j) = round(max(column) - min(column), 4);
    fprintf(num2str(differences(j)));
    fprintf(' ');
end
fprintf('\n');
optimal_degree = 1; % По умолчанию, если массив differences не содержит убывающей последовательности, то оптимальная степень - 1
for i = 1:length(differences)
    if differences(i) < differences(i+1) || differences(i) == 0.0
        break;
    else
        optimal_degree = i;
    end
end
disp(['Оптимальная степень полинома:', num2str(optimal_degree)]);
disp(['Ошибка в оценочной точке между 14 и 15 узлами: ', num2str(error1)]);
disp(['Ошибка в оценочной точке между 1 и 2 узлами: ', num2str(error2)]);

% Построение графиков
figure; 
plot(pressure, mass_heatCapacity, 'o', pressure, interpolation , '-', interp_temp, s, '--'); 
hold on;
plot(pressure, pl, '-'); 
hold on; 
plot(midpoint1, interp_value1, 'xr', midpoint2, interp_value2, 'xr'); 
xlabel('Давление (кПа)'); 
ylabel('Массовая теплоемкость'); 
legend('Узловые точки', 'Полином polyfit', 'Сплайн', 'Канонический полином' ,'Оценочные точки');
function [p] = buildCanonicalPoly(x, y)
    n = length(x);
    A = zeros(n, n);
    for i = 1:n
        A(:, i) = x.^(n-i);
    end
    p = det(vander(x)) * A\y';
end