% задание 4 варинт 14
% оперделяем матрицу соединений 
% Na N C H O Ca
A = [2 0 1 0 3 0; 
    0 1 0 1 3 0; 
    1 1 0 0 3 0; 
    0 0 0 2 1 0;
    0 0 1 0 2 0;
    0 0 0 0 1 1; 
    0 2 0 0 6 1];

disp(['Ранг матрицы: ', num2str(rank(A))]);
disp(['норма матрицы: ', num2str(norm(A, inf))]);
disp(['число обусловленности матрицы: ', num2str(cond(A))]);

% Инициализируем ячейку для хранения подматриц
submatrices = cell(1, (size(A, 1) - submatrixSize + 1) * (size(A, 2) - submatrixSize + 1));

% Извлекаем все подматрицы размером(rank*rank)
index = 1;
for i = 1:size(A, 1) - submatrixSize + 1
    for j = 1:size(A, 2) - submatrixSize + 1
        submatrices{index} = A(i:i+submatrixSize-1, j:j+submatrixSize-1);
        index = index + 1;
    end
end

% Печать всех подматриц
for k = 1:length(submatrices)
    det_A = det(submatrices{k});
    if abs(det_A) > 10^(-3)
        disp('невырожденный определитель != 10^(-3)');
         disp(submatrices{k});
    else
        disp('вырожденный определитель = 0 ');
         disp(submatrices{k});
    end
    disp('-----------');
end

%Количество возможных реакций
colSol = matrixSize(1)-submatrixSize;
disp(['Количество возможных реакций: ', num2str(colSol)]);

% Создание пустого массива Количество возможных реакций x размер столбца
B = repmat(0, colSol, matrixSize(1));

% Заполнение первых (Количество возможных реакций) элементов единичной матрицей
for i = 1:colSol
    B(i, i) = 1;
end
% Заполнение оставшихся элементов X с номером в матрице
for i = 1:colSol
    for j = (colSol + 1):matrixSize(1)
            B(i, j) = 1;
    end
end
disp('Матрица B');
disp(B);
disp('Матрица A');
disp(A);

% Инициализируем двумерный массив для хранения результатов умножения
result = cell(size(B, 1), size(A, 2));

% Перемножаем матрицы и записываем результат в массив
for i = 1:size(B, 1)
    for j = 1:size(A, 2)
        result_ij = zeros(size(B, 2), 1); % Инициализируем массив для хранения результата умножения
        for k = 1:size(B, 2)
            result_ij(k) =(B(i, k) * A(k, j)); % Выполняем побитовое "И" для каждой пары элементов
        end
        result{i, j} = result_ij; % Сохраняем результат умножения в соответствующий элемент массива
    end
end

for i = 1:colSol
    m = [transpose(result{i, 1}(end-4:end))];
    sm = [-sum(transpose(result{i, 1}(1:3)))];
    for j = 2:submatrixSize
        m = [m;transpose(result{1, j}(end-4:end))];
        sm = [sm;-sum(transpose(result{1, j}(1:3)))];
    end
    if i == 1
        sm_resultat = [sm];
        resultat = [m];
        X = round(linsolve(resultat, sm_resultat)); 
        error = resultat * X - sm_resultat; 
        disp(['ранг матрицы коэффициентов: ', num2str(rank(resultat))]);
        disp(['норма матрицы коэффициентов: ', num2str(norm(resultat, inf))]);
        disp(['число обусловленности матрицы коэффициентов: ', num2str(cond(resultat))]);
        disp('Решение:');
        disp(X);
        disp('Ошибка решения:');
        disp(error);

    else
        resultat = [resultat;m];
        sm_resultat = [sm_resultat;sm];
        X = round(linsolve(resultat, sm_resultat));
        disp(['ранг матрицы коэффициентов: ', num2str(rank(resultat))]);
        disp(['норма матрицы коэффициентов: ', num2str(norm(resultat, inf))]);
        disp(['число обусловленности матрицы коэффициентов: ', num2str(cond(resultat))]);
        disp('Решение:');
        disp(X);
        disp('Ошибка решения:');
        disp(error);
    end
end