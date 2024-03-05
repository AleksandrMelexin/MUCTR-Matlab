%1 
% ���� ������� ������������� A
A = [6, -1, -1;
     1, -2, 3;
     3, 4, 4];

% ���� ������� ��������� ������ B
B = [0; 1; -1];

% ������� �������� ������� invA
invA = inv(A);

% ������� ������� X
X = invA * B;

% ������� ���������
disp('resh:');
disp(X);

% ��������� ����������� ������� A
detA = det(A);

% ��������� ���� ������� A
rankA = rank(A);

% ��������� ����� ������� A
normA = norm(A);

% ��������� ����� ��������������� ������� A
condA = cond(A);

% ������� ����������
disp('det A:');
disp(detA);

disp('rank A:');
disp(rankA);

disp('norm A:');
disp(normA);

disp('chislo A:');
disp(condA);

disp('-------------------------------------------------------');
%--------------------------------------------------------------
%2
A = [9.1, 5.6, 7.8;3.8, 5.1, 2.8; 4.1, 5.7, 1.2];
b = [-9.8; 6.7; 5.8];
augmentedMatrix = [A, b];

rrefMatrix = rref(augmentedMatrix);
solution = rrefMatrix(:, end);

disp(solution);
% ��� ������������, ����� � ����� ��������������� �������� ����� �� ����������, ��� � ���������� �������

% ��������� ����������� ������� A
detA = det(A);

% ��������� ���� ������� A
rankA = rank(A);

% ��������� ����� ������� A
normA = norm(A);

% ��������� ����� ��������������� ������� A
condA = cond(A);

% ������� ����������
disp('det A:');
disp(detA);

disp('rank A:');
disp(rankA);

disp('norm A:');
disp(normA);

disp('chislo A:');
disp(condA);

disp('-------------------------------------------------------');
%3
% ���� ������� ������������� A
A = [2.34, -1.42, -0.54, 0.21;
     1.44, -0.53, 1.43, -1.27;
     0.63, -1.32, -0.65, 1.43;
     0.56, 0.88, -0.67, -2.38];

% ���� ������� ��������� ������ B
B = [0.66; -1.44; 0.94; 0.73];

% ������ ������� ��������� � ������� LU-����������
[L, U, P] = lu(A);
Y = L \ (P * B);
X = U \ Y;

% ������� ���������
disp('resh:');
disp(X);

% ��������� ����������� ������� A
detA = det(A);

% ��������� ���� ������� A
rankA = rank(A);

% ��������� ����� ������� A
normA = norm(A);

% ��������� ����� ��������������� ������� A
condA = cond(A);

% ������� ����������
disp('det A:');
disp(detA);

disp('rank A:');
disp(rankA);

disp('norm A:');
disp(normA);

disp('chislo A:');
disp(condA);

disp('-------------------------------------------------------');

%������� � 19
%����������� ������� ������� �� ��������� ������������ �����:
%H3PO4;  NH3;  NH4H2PO4;  (NH4)2HPO4;  H2SO4;  (NH4)2SO4;  H2O;   NH4OH.
Sostav = 'P  N  H  O  S'
P = [1 0 1 1 0 0 0 0];
N = [0 1 1 2 0 2 0 1];
H = [3 3 6 9 2 8 2 5];
O = [4 0 4 4 4 4 1 1];
S = [0 0 0 0 1 2 0 0];

% ����������� �������, ���������� ������ ������� �������
matrix = [P;N;H;O;S];
A = transpose(matrix)
matrixSize = size(A)
disp('-----------');

% ���������� ����� �������
submatrixSize = rank(A)

% �������������� ������ ��� �������� ���������
submatrices = cell(1, (size(A, 1) - submatrixSize + 1) * (size(A, 2) - submatrixSize + 1));

% ��������� ��� ���������� ��������(rank*rank)
index = 1;
for i = 1:size(A, 1) - submatrixSize + 1
    for j = 1:size(A, 2) - submatrixSize + 1
        submatrices{index} = A(i:i+submatrixSize-1, j:j+submatrixSize-1);
        index = index + 1;
    end
end

% ������ ���� ���������
for k = 1:length(submatrices)
    det_A = det(submatrices{k});
    if abs(det_A) > 10^(-3)
        disp('non - degenerate det <> 10^(-3)');
         disp(submatrices{k});
    else
        disp('degenerate det = 0 ');
         disp(submatrices{k});
    end
    disp('-----------');
end

%���������� ��������� �������
%num2str(matrixSize(2)
colSol = matrixSize(1)-submatrixSize;
disp(['The number of possible reactions ', num2str(colSol)]);

% �������� ������� ������� ���������� ��������� ������� x ������ �������
B = repmat(0, colSol, matrixSize(1));

% ���������� ������ (���������� ��������� �������) ��������� ��������� ��������
for i = 1:colSol
    B(i, i) = 1;
end
% ���������� ���������� ��������� X � ������� � �������
for i = 1:colSol
    for j = (colSol + 1):matrixSize(1)
            B(i, j) = 1;
    end
end
disp('the matrix B');
disp(B);
disp('the matrix A');
disp(A);

% �������������� ��������� ������ ��� �������� ����������� ���������
result = cell(size(B, 1), size(A, 2));

% ����������� ������� � ���������� ��������� � ������
for i = 1:size(B, 1)
    for j = 1:size(A, 2)
        result_ij = zeros(size(B, 2), 1); % �������������� ������ ��� �������� ���������� ���������
        for k = 1:size(B, 2)
            result_ij(k) =(B(i, k) * A(k, j)); % ��������� ��������� "�" ��� ������ ���� ���������
        end
        result{i, j} = result_ij; % ��������� ��������� ��������� � ��������������� ������� �������
    end
end


% ������� ���������
disp(result);


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
        % X = [linsolve(resultat, sm_resultat)]
         round(linsolve(resultat, sm_resultat))

    else
        resultat = [resultat;m];
        sm_resultat = [sm_resultat;sm];
        % X = [X;linsolve(resultat, sm_resultat)]
         round(linsolve(resultat, sm_resultat))
    end
end


