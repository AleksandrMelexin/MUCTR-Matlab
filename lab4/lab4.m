clear; clc;
A = [-0.76, -0.04, 0.21, -0.18; 0.45, -1.23, 0.66, 0; 0.26, 0.34, -1.11, 0; 0.05, -0.26, 0.34, -1.12]; % ��� �������� ����������������� �������
B = [-1.24; 0.88; -0.63; 1.17]; % ��� �������� ����������������� �������
% ������� 14
%A = [8 4 -6 0; 1 2 1 -6; -3 -6 -2 -9; 4 3 2 1];
%B = [596; 262.02; -731.47; 396.83];

% 1 ����������� ������������ ������� �������������
det_A = det(A);
disp(['����������� ������� �������������: ', num2str(det_A)]);

% 2 ���� ������� �������������, �����, ����� ���������������
rank_A = rank(A);
norm_A = norm(A);
cond_A = cond(A);

disp(['���� ������� �������������: ', num2str(rank_A)]);
disp(['����� ������� �������������: ', num2str(norm_A)]);
disp(['����� ��������������� ������� �������������: ', num2str(cond_A)]);

fprintf('\n---------------------------------------------------------\n\n');
% 3 ������� �������� ������� �������
eps = 1e-3;
disp(['�������� ������� �������: ', num2str(eps)]);

fprintf('\n---------------------------------------------------------\n\n');

% 4 ������ ������� ������� ������� ��������
% k - ������������ ���������� �������� (����� ����������� ���� ����� ����)
k = 1000;
% �������� ������� ���������� ��� ������ ������� ��������
T_simple = inv(diag(diag(A))) * (diag(diag(A)) - A);
spectral_radius_simple = max(abs(eig(T_simple)));
% ������ ��� - �������� ���������� �������� � ����������� ��������
if spectral_radius_simple >= 1
    disp('����� ������� �������� ����������');
else
    disp('����� ������� �������� ��������');
    t = 0.5; % ��������� ��� ������ ������� ��������
    [X_simple, iters_simple] = simple_iteration_method(A, B, k, eps, t); 
    disp('������� ������� ������� ��������:');
    disp(X_simple);
    disp(['����� �������� ��� �������: ', num2str(iters_simple)]);
end

fprintf('\n---------------------------------------------------------\n\n');

% 5. ������� ������� ������� �������

% �������� ������� ���������� ��� ������ �������
k = 1000;
[X_seidel, iters_seidel] = seidel_method(A, B, k, eps);
if ~any(~isnan(X_seidel(:)))
    disp('����� ������� ����������');
else
    disp('����� ������� ��������');
    disp('������� ������� �������:');
    disp(X_seidel);
    disp(['����� ��������� ��������: ', num2str(iters_seidel)]);
end

fprintf('\n---------------------------------------------------------\n\n');

% 6 ������� ������� ������� �����
% �������� ������� ���������� ��� ������ �����
%if isdiagonaldominant(A)
[X_jacobi, iters_jacobi] = jacobi_method(A, B, k, eps);
if ~any(~isnan(X_jacobi(:)))
    disp('����� ����� ����������');
else
    disp('����� ����� ��������');
    disp('������� ������� �����:');
    disp(X_jacobi);
    disp(['����� ��������� ��������: ', num2str(iters_jacobi)]);
end    

fprintf('\n---------------------------------------------------------\n\n');

disp('������� linsolve()')
X = linsolve(A, B);
disp(transpose(X));

disp('�������: ')
for i = 1:rank_A
    val = sprintf('x%u = %f \n', i, X(i));
    fprintf(val);
end

fprintf('\n---------------------------------------------------------\n\n');

function [X, iters] = seidel_method(A, B, max_iters, epsilon) % ����� �������
    X = zeros(size(B)); % ��������� ����������� 
    for iters = 1:max_iters
        X_old = X;
        for i = 1:length(B)
            sigma = A(i, 1:i-1) * X(1:i-1) + A(i, i+1:end) * X_old(i+1:end);
            X(i) = (B(i) - sigma) / A(i, i);
        end
        if norm(X - X_old, inf) < epsilon % �������� ������� ���������(����� �������� ������ ��������)
            break;
        end
    end
end

function [X, iters] = jacobi_method(A, B, max_iters, epsilon) % ����� �����
    X = zeros(size(B)); % ��������� ����������� 
    n = length(B);
    for iters = 1:max_iters
        X_old = X;
        for i = 1:n
            % ����� ���� ��������� ������ ������� A, 
            % ���������� �� ��������������� �������� X_old, 
            % �� ����������� ������������� ��������, ������� ���������� �� �����, ����������� �� ��������������� �������� X_old.
            sigma = A(i, :) * X_old - A(i, i) * X_old(i);
            X(i) = (B(i) - sigma) / A(i, i);
        end
        if norm(X - X_old, inf) < epsilon % �������� ������� ���������(����� �������� ������ ��������)
            break;
        end
    end
end

function [X, iters] = simple_iteration_method(A, B, max_iters, epsilon, t) % ����� ������� ��������
    X = zeros(size(B)); % ��������� �����������
    % max_iters - ������������ ���������� �������� (����� ����������� ���� ����� ����)
    for iters = 1:max_iters
        X_new = A * X + t*B; % ���������������� ���������� ������ �����������
        if norm(X_new - X, inf) < epsilon % �������� ������� ���������(����� �������� ������ ��������)
            X = X_new;
            break;
        end
        X = X_new;
    end
end


