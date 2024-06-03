clc; clear;

% Параметры
tol = 1e-2; % Увеличение допустимой погрешности для уменьшения частоты уменьшения шага
h = 0.1;

% Инициализация для явного метода Эйлера
x_explicit = 0;
y1_explicit = 1;
y2_explicit = 1;

% Массивы для хранения чисел жёсткости и шагов
stiffness_numbers = [];
steps = [];

i = 1;
while x_explicit(end) <= 1
    % Матрица Якоби и число жёсткости
    g = [exp(x_explicit(end)^2), x_explicit(end); -1, 2];
    s = max(real(eig(g))) / min(real(eig(g)));
    stiffness_numbers(end + 1) = s;

    % Оценка ошибок для явного метода Эйлера
    y1_temp = y1_explicit(end) + (y1_explicit(end) * exp(-x_explicit(end)^2) + x_explicit(end) * y2_explicit(end)) * h;
    y2_temp = y2_explicit(end) + (3 * x_explicit(end) - y1_explicit(end) + 2 * y2_explicit(end)) * h;

    y1_explicit_half = y1_explicit(end) + (y1_explicit(end) * exp(-x_explicit(end)^2) + x_explicit(end) * y2_explicit(end)) * (h / 2);
    y2_explicit_half = y2_explicit(end) + (3 * x_explicit(end) - y1_explicit(end) + 2 * y2_explicit(end)) * (h / 2);

    y1_explicit_double = y1_explicit_half + (y1_explicit_half * exp(-(x_explicit(end) + h / 2)^2) + (x_explicit(end) + h / 2) * y2_explicit_half) * (h / 2);
    y2_explicit_double = y2_explicit_half + (3 * (x_explicit(end) + h / 2) - y1_explicit_half + 2 * y2_explicit_half) * (h / 2);

    error1 = abs(y1_explicit_double - y1_temp);
    error2 = abs(y2_explicit_double - y2_temp);
    error = max(error1, error2);

    if error > tol
        h = h / 2;  % уменьшение шага
    elseif error < tol / 4
        h = h * 2;  % увеличение шага
    end

    % Принятие шага
    y1_explicit(end + 1) = y1_temp;
    y2_explicit(end + 1) = y2_temp;
    x_explicit(end + 1) = x_explicit(end) + h;
    
    % Запись текущего шага в массив
    steps(end + 1) = h;
    i = i + 1;
end

% Вывод чисел жёсткости для явного йлера
fprintf('Числа жёсткости:\n');
for i = 1:length(stiffness_numbers)
    fprintf('Итерация %d: число жёсткости = %f\n', i, stiffness_numbers(i));
end

fprintf('\n\n');
% Вывод шагов для явного йлера
fprintf('Шаги:\n');
for i = 1:length(steps)
    fprintf('Итерация %d: шаг = %f\n', i, steps(i));
end

fprintf('\n\n');
% Вывод найденных значений функции
fprintf('Значения функции y1 = f1(x) по явному методу Эйлера:\n\n');
for i = 1:length(x_explicit)
    disp([num2str(x_explicit(i)), '   ', num2str(y1_explicit(i))]);
end
fprintf('\n');
fprintf('Значения функции y2 = f2(x) по явному методу Эйлера:\n\n');
for i = 1:length(x_explicit)
    disp([num2str(x_explicit(i)), '   ', num2str(y2_explicit(i))]);
end
disp('___________________________________________________________________')

% Инициализация для неявного метода Эйлера
x_implicit = 0;
y1_implicit = 1;
y2_implicit = 1;

i = 1;

while x_implicit(end) <= 1
    % Матрица Якоби и число жёсткости
    g = [exp(x_explicit(end)^2), x_explicit(end); -1, 2];
    s = max(real(eig(g))) / min(real(eig(g)));
    stiffness_numbers(end + 1) = s;

    % Оценка ошибок для неявного метода Эйлера
    y1_temp = y1_implicit(end) + (y1_implicit(end) * exp(x_implicit(end)^2) + x_implicit(end) * y2_implicit(end)) * h;
    y2_temp = y2_implicit(end) + (3 * x_implicit(end) - y1_implicit(end) + 2 * y2_implicit(end)) * h;

    y1_half = y1_implicit(end) + (y1_implicit(end) * exp(x_implicit(end)^2) + x_implicit(end) * y2_implicit(end) + y1_temp) * (h / 2);
    y2_half = y2_implicit(end) + (3 * x_implicit(end) - y1_implicit(end) + 2 * y2_implicit(end) + 3 * x_implicit(end) - y1_temp + 2 * y2_temp) * (h / 2);

    error1 = abs(y1_half - y1_temp);
    error2 = abs(y2_half - y2_temp);
    error = max(error1, error2);

    if error > tol
        h = h / 2;  % уменьшение шага
    elseif error < tol / 4
        h = h * 2;  % увеличение шага
    end

    % Принятие шага
    y1_implicit(end + 1) = y1_temp;
    y2_implicit(end + 1) = y2_temp;
    x_implicit(end + 1) = x_implicit(end) + h;
    
    % Запись текущего шага в массив
    steps(end + 1) = h;
    i = i + 1;
end

% Вывод найденных значений функции
fprintf('Значения функции y1 = f1(x) по неявному методу Эйлера:\n\n');
for i = 1:length(x_implicit)
    disp([num2str(x_implicit(i)), '   ', num2str(y1_implicit(i))]);
end
fprintf('Значения функции y2 = f2(x) по неявному методу Эйлера:\n\n');
for i = 1:length(x_implicit)
    disp([num2str(x_implicit(i)), '   ', num2str(y2_implicit(i))]);
end
disp('___________________________________________________________________')

% Стандартные методы ODE15s
% Определение начальных условий
y0 = [1; 1];

% Решение системы дифференциальных уравнений
[x_ode15s, y_ode15s] = ode15s(@myODEs, [0, 1], y0);

% Вывод чисел жёсткости для неявного Эйлера
fprintf('Числа жёсткости:\n');
for i = 1:length(stiffness_numbers)
    fprintf('Итерация %d: число жёсткости = %f\n', i, stiffness_numbers(i));
end

% Вывод шагов для неявного йлера
fprintf('Шаги:\n');
for i = 1:length(steps)
    fprintf('Итерация %d: шаг = %f\n', i, steps(i));
end

% Визуализация результатов
figure
hold on
plot(x_explicit, y1_explicit, 'r', 'DisplayName', 'Явный Эйлер y1');
plot(x_explicit, y2_explicit, 'g', 'DisplayName', 'Явный Эйлер y2');
plot(x_implicit, y1_implicit, 'b', 'DisplayName', 'Неявный Эйлер y1');
plot(x_implicit, y2_implicit, 'm', 'DisplayName', 'Неявный Эйлер y2');
plot(x_ode15s, y_ode15s(:, 1), 'k--', 'DisplayName', 'MATLAB y1 (ODE15s)');
plot(x_ode15s, y_ode15s(:, 2), 'c--', 'DisplayName', 'MATLAB y2 (ODE15s)');
grid on
xlabel('x');
ylabel('y');
legend;
title('Решение системы дифференциальных уравнений');
hold off

% Определение функции myODEs
function dydx = myODEs(x, y)
    dydx = [y(1) * exp(x^2) + x * y(2); 3 * x - y(1) + 2 * y(2)];
end