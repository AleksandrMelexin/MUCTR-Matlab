clear; clc;
A = [-0.76, -0.04, 0.21, -0.18; 0.45, -1.23, 0.66, 0; 0.26, 0.34, -1.11, 0; 0.05, -0.26, 0.34, -1.12]; % для проверки работоспособности методов
B = [-1.24; 0.88; -0.63; 1.17]; % для проверки работоспособности методов
% вариант 14
%A = [8 4 -6 0; 1 2 1 -6; -3 -6 -2 -9; 4 3 2 1];
%B = [596; 262.02; -731.47; 396.83];

% 1 Определение детерминанта матрицы коэффициентов
det_A = det(A);
disp(['Детерминант матрицы коэффициентов: ', num2str(det_A)]);

% 2 Ранг матрицы коэффициентов, норма, число обусловленности
rank_A = rank(A);
norm_A = norm(A);
cond_A = cond(A);

disp(['Ранг матрицы коэффициентов: ', num2str(rank_A)]);
disp(['Норма матрицы коэффициентов: ', num2str(norm_A)]);
disp(['Число обусловленности матрицы коэффициентов: ', num2str(cond_A)]);

fprintf('\n---------------------------------------------------------\n\n');
% 3 Задание точности решения системы
eps = 1e-3;
disp(['Точность решения системы: ', num2str(eps)]);

fprintf('\n---------------------------------------------------------\n\n');

% 4 Решить систему методом простых итераций
% k - максимальное количество итераций (иначе бесконечный цикл может быть)
k = 1000;
% Проверка условий сходимости для метода простой итерации
T_simple = inv(diag(diag(A))) * (diag(diag(A)) - A);
spectral_radius_simple = max(abs(eig(T_simple)));
% спектр рад - максимум абсолютных значений её собственных значений
if spectral_radius_simple >= 1
    disp('Метод простых итераций расходится');
else
    disp('Метод простых итераций сходится');
    t = 0.5; % множитель для метода простых итераций
    [X_simple, iters_simple] = simple_iteration_method(A, B, k, eps, t); 
    disp('Решение методом простых итераций:');
    disp(X_simple);
    disp(['Число итераций для решения: ', num2str(iters_simple)]);
end

fprintf('\n---------------------------------------------------------\n\n');

% 5. Решение системы методом Зейделя

% Проверка условий сходимости для метода Зейделя
k = 1000;
[X_seidel, iters_seidel] = seidel_method(A, B, k, eps);
if ~any(~isnan(X_seidel(:)))
    disp('Метод Зейделя расходится');
else
    disp('Метод Зейделя сходится');
    disp('Решение методом Зейделя:');
    disp(X_seidel);
    disp(['Число требуемых итераций: ', num2str(iters_seidel)]);
end

fprintf('\n---------------------------------------------------------\n\n');

% 6 Решение систему методом Якоби
% Проверка условий сходимости для метода Якоби
%if isdiagonaldominant(A)
[X_jacobi, iters_jacobi] = jacobi_method(A, B, k, eps);
if ~any(~isnan(X_jacobi(:)))
    disp('Метод Якоби расходится');
else
    disp('Метод Якоби сходится');
    disp('Решение методом Якоби:');
    disp(X_jacobi);
    disp(['Число требуемых итераций: ', num2str(iters_jacobi)]);
end    

fprintf('\n---------------------------------------------------------\n\n');

disp('Решение linsolve()')
X = linsolve(A, B);
disp(transpose(X));

disp('Решение: ')
for i = 1:rank_A
    val = sprintf('x%u = %f \n', i, X(i));
    fprintf(val);
end

fprintf('\n---------------------------------------------------------\n\n');

function [X, iters] = seidel_method(A, B, max_iters, epsilon) % метод зейделя
    X = zeros(size(B)); % Начальное приближение 
    for iters = 1:max_iters
        X_old = X;
        for i = 1:length(B)
            sigma = A(i, 1:i-1) * X(1:i-1) + A(i, i+1:end) * X_old(i+1:end);
            X(i) = (B(i) - sigma) / A(i, i);
        end
        if norm(X - X_old, inf) < epsilon % проверка условия остановки(норма разности меньше точности)
            break;
        end
    end
end

function [X, iters] = jacobi_method(A, B, max_iters, epsilon) % метод якоби
    X = zeros(size(B)); % Начальное приближение 
    n = length(B);
    for iters = 1:max_iters
        X_old = X;
        for i = 1:n
            % сумма всех элементов строки матрицы A, 
            % умноженных на соответствующие значения X_old, 
            % за исключением диагонального элемента, который вычитается из суммы, домноженный на соответствующее значение X_old.
            sigma = A(i, :) * X_old - A(i, i) * X_old(i);
            X(i) = (B(i) - sigma) / A(i, i);
        end
        if norm(X - X_old, inf) < epsilon % проверка условия остановки(норма разности меньше точности)
            break;
        end
    end
end

function [X, iters] = simple_iteration_method(A, B, max_iters, epsilon, t) % метод простых итераций
    X = zeros(size(B)); % Начальное приближение
    % max_iters - максимальное количество итераций (иначе бесконечный цикл может быть)
    for iters = 1:max_iters
        X_new = A * X + t*B; % последовательное нахождение нового приближения
        if norm(X_new - X, inf) < epsilon % проверка условия остановки(норма разности меньше точности)
            X = X_new;
            break;
        end
        X = X_new;
    end
end


