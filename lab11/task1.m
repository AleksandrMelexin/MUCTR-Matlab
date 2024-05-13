clc; clear;
% Задание 1
% Начальный интервал
a = 10;
b = 30;
d0 = 5; % Начальное приближение
tol = 1e-4; % Точность и максимальное количество итераций
max_iter = 100;

disp('----------------------------------------------------------------------------------')
disp('Задание 1')

% Вызов методов Ньютона, Золотого чесения и парабол
min_d_newton = newton_method(d0, tol, max_iter);
min_d_golden = golden_section_method(a, b, tol, max_iter);
min_d_parabolic = parabolic_method(a, b, tol, max_iter);

% Вызов функции MATLAB
options = optimset('Display', 'off'); 
[min_d_matlab, min_heat_loss] = fminsearch(@calculate_heat_loss, d0, options);

% Вывод результата
disp(['Оптимальное значение d (метод Ньютона): ', num2str(min_d_newton)]);
disp(['Оптимальное значение d (метод золотого сечения): ', num2str(min_d_golden)]);
disp(['Оптимальное значение d (метод парабол): ', num2str(min_d_parabolic)]);
disp(['Оптимальное значение d (стандартная функция MATLAB): ', num2str(min_d_matlab)]);
disp('----------------------------------------------------------------------------------')

function heat_loss = calculate_heat_loss(d) % Функция для вычисления потерь тепла
    h = (4*1200)./(pi*d.^2); % выражение высоты исходя из общего объёма
    % Коэффициенты теплоотдачи
    h_bottom = 0.35; % Для дна
    h_wall = 1; % Для стенок
    h_lid = 0.68; % Для крышки
    
    % Площади поверхностей
    S_bottom = pi * (d.^2) / 4;
    S_wall = pi * d .* h;
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
    iter = 0; % Инициализация счетчика итераций
    d_values = []; % Для хранения значений d на каждой итерации
    heat_losses = []; % Для хранения значений потерь тепла на каждой итерации
    while iter < max_iter
        h_loss = calculate_heat_loss(d);
        h = 1e-4; % малое приращение для численного дифференцирования
        d_prime = (calculate_heat_loss(d + h) - calculate_heat_loss(d)) / h;
        d_double_prime = (calculate_heat_loss(d + h) - 2 * h_loss + calculate_heat_loss(d - h)) / h^2;
        d = d - d_prime / d_double_prime;
        % Сохраняем значения d и потерь тепла на текущей итерации
        d_values = [d_values, d];
        heat_losses = [heat_losses, h_loss];
        if abs(d_prime) < tol
            break;
        end
        iter = iter + 1; % Увеличение счетчика итераций
    end
    min_d = d;
    d_values(end) = [];
    heat_losses(end) = [];
    disp(['Количество итераций (метод Ньютона): ', num2str(iter)]); % Вывод количества итераций на экран

    % Вывод графика
    figure;
    plot(linspace(10, 30, 100), calculate_heat_loss(linspace(10, 30, 100)), 'b'); % График функции потерь тепла
    hold on;
    plot(min_d, calculate_heat_loss(min_d), 'ro'); % Точка экстремума
    plot(d_values, heat_losses, 'ko'); % Текущие значения d
    xlabel('d');
    ylabel('Потери тепла');
    title('Метод Ньютона');
    legend('Потери тепла', 'Экстремум', 'Текущие значения d');
    grid on;
    hold off;
end

function min_d = golden_section_method(a, b, tol, max_iter) % метод золотого сечения
    rho = (sqrt(5) - 1) / 2;
    d = a + rho * (b - a);
    c = b - rho * (b - a);
    iter = 0; % Инициализация переменной iter
    d_values = []; % Для хранения значений d на каждой итерации
    heat_losses = []; % Для хранения значений потерь тепла на каждой итерации
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
        % Сохраняем значения d и потерь тепла на текущей итерации
        d_values = [d_values, d];
        heat_losses = [heat_losses, calculate_heat_loss(d)];
        iter = iter + 1; % Увеличение счетчика итераций
    end
    min_d = (a + b) / 2;
    d_values(end) = [];
    heat_losses(end) = [];
    disp(['Количество итераций (метод золотого сечения): ', num2str(iter)]); % Вывод количества итераций на экран

    % Вывод графика
    figure;
    plot(linspace(10, 30, 100), calculate_heat_loss(linspace(10, 30, 100)), 'b'); % График функции потерь тепла
    hold on;
    plot(min_d, calculate_heat_loss(min_d), 'ro'); % Точка экстремума
    plot(d_values, heat_losses, 'ko'); % Текущие значения d
    xlabel('d');
    ylabel('Потери тепла');
    title('Метод золотого сечения');
    legend('Потери тепла', 'Экстремум', 'Текущие значения d');
    grid on;
    hold off;
end

function min_d = parabolic_method(a, b, tol, max_iter) % метод парабол
    x = [a, (a + b) / 2, b];
    f = zeros(1, 3);
    f(1) = calculate_heat_loss(x(1));
    f(2) = calculate_heat_loss(x(2));
    f(3) = calculate_heat_loss(x(3));
    iter = 0;
    d_values = []; % Для хранения значений d на каждой итерации
    heat_losses = []; % Для хранения значений потерь тепла на каждой итерации
    while abs(x(2) - x(1)) > tol && iter < max_iter
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
        % Сохраняем значения d и потерь тепла на текущей итерации
        d_values = [d_values, x(2)];
        heat_losses = [heat_losses, f(2)];
        iter = iter + 1;
    end
    min_d = x(2);
    d_values(end) = [];
    heat_losses(end) = [];
    disp(['Количество итераций (метод парабол): ', num2str(iter)]); % Вывод количества итераций на экран

    % Вывод графика
    figure;
    plot(linspace(10, 30, 100), calculate_heat_loss(linspace(10, 30, 100)), 'b'); % График функции потерь тепла
    hold on;
    plot(min_d, calculate_heat_loss(min_d), 'ro'); % Точка экстремума
    plot(d_values, heat_losses, 'ko'); % Текущие значения d
    xlabel('d');
    ylabel('Потери тепла');
    title('Метод парабол');
    legend('Потери тепла', 'Экстремум', 'Текущие значения d');
    grid on;
    hold off;
end
