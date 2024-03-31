clear; clc; % ������� ���������� ���� � ����������
% ����������� ������� ��� ��������� ���������
f = @(x) 2 * exp(-3 * x) - 3 * x + 2;

% ����������� ��������� ��� �������� ������
a = -2; % ������ ���������
b = 2; % ����� ���������
n = 1000; % ���������� ����� ��� ��������

% ����� ��������
x = linspace(a, b, n); % ��������� n ���������� �������������� ����� �� ��������� [a, b]
y = f(x); % ���������� �������� ������� � ���� ������
roots = x(abs(y) < 0.01); % ����� ������ � ������������ ���������

% ���������� ��������
disp(' ');
fprintf('-------------------------\n\n');
fprintf('����� ��������� (����� ��������): ');
disp(roots);
fprintf('-------------------------\n\n');

% ������ ������� � ��������� ������
figure;
plot(x, y, 'b'); % ���������� ������� �������
hold on;
plot(roots, f(roots), 'm*'); % ������� ������ �� �������
xlabel('x');
ylabel('f(x)');
title('������ ������� � ��������� ������');
legend('������� f(x)', '�����', 'Location', 'NorthWest');

% �������� ���������� ������ ����������� �������

% �������� ���������� ����������� a � b
a = -1;
b = 1;
tol = 1e-6;
maxIter = 100;

% �������� ������� ����������
converge = true; % ������������, ��� ��������
iter = 0;
if f(a) * f(b) > 0
    converge = false; % �� ��������
end

while abs(f(a) * f(b)) > tol && iter < maxIter
    c = (a + b) / 2;
    if f(a) * f(c) < 0
        b = c;
    else
        a = c;
    end
    iter = iter + 1;
end

fprintf('���������� ������ ����������� �������: %d (1 - ��������, 0 - �� ��������)\n', converge);


% ����� ����������� �������
a = -1; 
b = 1; 
tol = 1e-6; % ��������
maxIter = 100; % ������������ ���������� ��������

if f(a) * f(b) > 0
    error('������� ���������� �� ���������. �������� ������ ��������.');
end

c = (a + b) / 2;  % ��������� �����������

iter = 0;
while abs(f(c)) > tol && iter < maxIter
    if f(a) * f(c) < 0
        b = c;
    else
        a = c;
    end
    c = (a + b) / 2;
    iter = iter + 1;
end
fprintf('����� ����������� �������: ������ %f, ����� ��������: %d\n', c, iter);
fprintf('-------------------------\n\n');
% �������� ���������� ������ ������� ��������
% ����������� ������� g(x)
dg = @(x) (3 / (3 * x - 2));

% �������� ������� ���������� |g'(x)| < 1 � ����������� �����
x = 1.8; % ����� ���������� �������� x � ����������� �����

maxIter = 100; % ������������ ���������� ��������

% �������� ������� ����������
converge = true; % ������������, ��� ��������
iter = 0;
while abs(dg(x)) >= 1 && iter < maxIter
    x = g(x);
    iter = iter + 1;
end
if abs(dg(x)) >= 1
    converge = false; % �� ��������
end

fprintf('���������� ������ ������� ��������: %d (1 - ��������, 0 - �� ��������)\n', converge);

% ����� ������� ��������
g = @(x) log((3 * x - 2) / 2);
x0 = 0; % ��������� �����������
tol = 1e-6; % ��������
maxIter = 100; % ������������ ���������� ��������

iter = 0;
while iter < maxIter
    x = g(x0);
    if abs(x - x0) < tol
        break;
    end
    x0 = x;
    iter = iter + 1;
end
fprintf('����� ������� ��������: ������ %f, ����� ��������: %d\n', x, iter);
fprintf('-------------------------\n\n');

% �������� ���������� ������ ����

% ��������� ����������� x0 � x1
x0 = -1;
x1 = 1;
tol = 1e-6; % ��������
maxIter = 100; % ������������ ���������� ��������

% �������� ������� ����������
converge = true; % ������������, ��� ��������
iter = 0;
while abs(x1 - x0) >= tol && iter < maxIter
    x = x1 - (f(x1) * (x1 - x0)) / (f(x1) - f(x0));
    if abs(x - x1) < tol
        break;
    end
    x0 = x1;
    x1 = x;
    iter = iter + 1;
end
if abs(x1 - x0) >= tol
    converge = false; % �� ��������
end

fprintf('���������� ������ ����: %d (1 - ��������, 0 - �� ��������)\n', converge);

% ����� ����
x0 = 0; % ��������� �����������
x1 = 1; % ������ ��������� �����������
tol = 1e-6; % ��������
maxIter = 100; % ������������ ���������� ��������

iter = 0;
while iter < maxIter
    x = x1 - (f(x1) * (x1 - x0)) / (f(x1) - f(x0));
    if abs(x - x1) < tol
        break;
    end
    x0 = x1;
    x1 = x;
    iter = iter + 1;
end
fprintf('����� ����: ������ %f, ����� ��������: %d\n', x, iter);
fprintf('-------------------------\n\n');

% �������� ���������� ������ �������

% ��������� ����������� x0 � x1
x0 = -1;
x1 = 1;
tol = 1e-6; % ��������
maxIter = 100; % ������������ ���������� ��������

% �������� ������� ����������
converge = true; % ������������, ��� ��������
iter = 0;
while abs(x1 - x0) >= tol && iter < maxIter
    x = x1 - (f(x1) * (x1 - x0)) / (f(x1) - f(x0));
    if abs(x - x1) < tol
        break;
    end
    x0 = x1;
    x1 = x;
    iter = iter + 1;
end
if abs(x1 - x0) >= tol
    converge = false; % �� ��������
end

fprintf('���������� ������ �������: %d (1 - ��������, 0 - �� ��������)\n', converge);





% ����� �������
x0 = -1; % ��������� ����������� 1
x1 = 1; % ��������� ����������� 2
tol = 1e-6; % ��������
maxIter = 100; % ������������ ���������� ��������

iter = 0;
while iter < maxIter
    x = x1 - (f(x1) * (x1 - x0)) / (f(x1) - f(x0));
    if abs(x - x1) < tol
        break;
    end
    x0 = x1;
    x1 = x;
    iter = iter + 1;
end
fprintf('����� �������: ������ %f, ����� ��������: %d\n', x, iter);
fprintf('-------------------------\n\n');


