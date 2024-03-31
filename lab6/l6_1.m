clear;
clc;

syms x y

eq1 = (x^2 + y^2 - 1)^3 - x^2 * y ^3 == 0;
eq2 = -y - 0.05 * x^2 - 0.5 * x + 0.25 == 0;

eqns = [eq1, eq2];
vars = [x, y];

f1 = @(x, y) (x^2 + y^2 - 1)^3 - x^2 * y ^3;
%fimplicit(f1, [-2, 2, -2, 2]);
ezplot(f1, [-2, 2]);
hold on
f2 = @(x, y) -y - 0.05 * x^2 - 0.5 * x + 0.25;
%fimplicit(f2, [-2, 2, -2, 2]);
ezplot(f2, [-2, 2]);
hold off

% Простые итерации
g1 = @(x, y) sqrt(1 - y^2);
g2 = @(x, y) x^2;
x0 = 0.8;
y0 = 0.6;
max_iter = 1000;
tol = 1e-6;

x = x0;
y = y0;
for i = 1:max_iter
    x_next = g1(x, y);
    y_next = g2(x, y);
    if abs(x_next - x) < tol && abs(y_next - y) < tol
        break;
    end
    x = x_next;
    y = y_next;
end

fprintf('Метод простых итераций:\n');
fprintf('x = %f, y = %f\n', x, y);
fprintf('Количество итераций: %d\n', i);


% Метод Ньютона
f = @(x, y) [x^2 + y^2 - 1; x^2 - y];
J = @(x, y) [2*x, 2*y; 2*x, -1];
x0 = 0.8;
y0 = 0.6;
max_iter = 100;
tol = 1e-6;



x = x0;
y = y0;
for i = 1:max_iter
    F = f(x, y);
    J_inv = inv(J(x, y));
    %disp(J_inv);
    %disp('---------------------------------------------');
    %disp(det(J_inv));
    d = J_inv * F;
    x = x - d(1);
    y = y - d(2);
    if norm(F) < tol
        break;
    end
end

fprintf('\nМетод Ньютона:\n');
fprintf('x = %f, y = %f\n', x, y);
fprintf('Количество итераций: %d\n', i);


%[solX, solY] = fsolve(@(vars) [vars(1)^2 + vars(2)^2 - 1; vars(1)^2 - vars(2)], [0.8, 0.6]);
[solX, solY] = fsolve(@(vars) [vars(1)^2 + vars(2)^2 - 1; vars(1)^2 - vars(2)], [0.8, 0.6]);
fprintf('\nС использованием функции fsolve:\n');
fprintf('x = %f, y = %f\n', solX, solY);

[solX, solY] = vpasolve(eqns, vars);
fprintf('\nС использованием vpasolve:\n');
disp(solX);
disp(solY);
fprintf('x = %f, y = %f\n', solX, solY);


%fimplicit(f1, [-2, 2, -2, 2]);
ezplot(f1, [-2, 2]);
hold on
%fimplicit(f2, [-2, 2, -2, 2]);
ezplot(f1, [-2, 2]);
plot(x, y, 'r*', 'MarkerSize', 10);
hold off