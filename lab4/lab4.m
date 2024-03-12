disp('вариант 14');
disp('задание 1');
A = [-0.76 0.04 0.21 0.18;
    0.45 -1.23 0.66 0;
    0.26 0.34 -1.11 0;
    0.05 -0.26 0.34 -1.12];
b = [-1.24;
    0.88;
    -0.63;
    1.17];
tol = 1e-6; % Погрешность
max_iter = 1000; % Максимальное число итераций

[x_simple, iter_simple, x_seidel, iter_seidel] = iterative_methods(A, b, tol, max_iter);

% Вывод результатов
disp('Метод простой итерации:');
disp(['Корни: ', num2str(x_simple')]);
disp(['Число итераций: ', num2str(iter_simple)]);

disp('Метод Зейделя:');
disp(['Корни: ', num2str(x_seidel')]);
disp(['Число итераций: ', num2str(iter_seidel)]);

% Проверка, подставив найденные корни в уравнения
disp('Проверка уравнений:');
disp('A * x_simple + b:');
disp(A * x_simple + b);
disp('A * x_seidel + b:');
disp(A * x_seidel + b);

% Решение уравнений с помощью стандартных операторов MATLAB
disp('Решение уравнений с помощью MATLAB:');
disp('x_sol:');
x_sol = A\b;
disp(x_sol);

%------------------------------------------------------------------
disp('задание 2');
A = [2 1 4 1;
    3 0 1 1;
    -1 2 3 4;
    3 1 1 1];

% Найдем собственные значения и собственные векторы
[V, D] = eig(A);

% V содержит собственные векторы в столбцах
% D содержит собственные значения на диагонали

% Выведем собственные значения и собственные векторы
disp('Собственные значения:');
disp(diag(D)');
disp('Собственные векторы:');
disp(V);

%------------------------------------------------------------------
disp('задание 3');
A = [8 4 -6 0;
    1 2 1 -6;
    -3 -6 -2 -9;
    4 3 2 1];
b = [596;
    262.02;
    -731.47;
    396.83];

function [x_simple, iter_simple, x_seidel, iter_seidel] = iterative_methods(A, b, tol, max_iter)
    f_i = 1;
    f_z = 1;
    n = length(b);
    x_simple = zeros(n, 1);
    x_seidel = zeros(n, 1);
    iter_simple = 0;
    iter_seidel = 0;
    
    % Проверка условия сходимости метода простых итераций
    if max(abs(eig(A))) >= 1
        disp('Метод простой итерации не сходится для данной матрицы');
    end
    
    % Проверка условия сходимости метода Зейделя
    if max(abs(eig(tril(A)))) >= 1
        disp('Метод Зейделя не сходится для данной матрицы');
    end
    
    % Метод простой итерации
    while true
        x_new = A * x_simple + b;
        if norm(x_new - x_simple, inf) < tol || iter_simple >= max_iter
            break;
        end
        x_simple = x_new;
        iter_simple = iter_simple + 1;
    end
    
    % Метод Зейделя
    while true
        x_new = zeros(n, 1);
        for i = 1:n
            x_new(i) = (b(i) - A(i, 1:i-1) * x_new(1:i-1) - A(i, i+1:end) * x_seidel(i+1:end)) / A(i, i);
        end
        if norm(x_new - x_seidel, inf) < tol || iter_seidel >= max_iter
            break;
        end
        x_seidel = x_new;
        iter_seidel = iter_seidel + 1;
    end
end