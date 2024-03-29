% задание 1
A = [6 -1 -1; % матрица коэффициентов
    1 -2 3;
    3 4 4];
B = [0; 1; -1]; % столбец свободных коэффициентов
X1 = inv(A)*B; % столбец неизвестных
X2 = linsolve(A, B);
error = A * X1 - B; % проверка точности решения
disp(['детерминант матрицы коэффициентов: ', num2str(det(A))]); 
disp(['ранг матрицы коэффициентов: ', num2str(rank(A))]);
disp(['норма матрицы коэффициентов: ', num2str(norm(A, inf))]);
disp(['число обусловленности матрицы коэффициентов: ', num2str(cond(A))]);
disp('Решение СЛАУ методом обратной матрицы:');
disp(X1);
disp('Ошибка решения:');
disp(error);
disp('Решение СЛАУ методом с помощью функции linsolve:');
disp(X2);

