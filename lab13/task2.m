clc; clear;

% Эйлер явный
x_explicit(1) = 0;
y1_explicit(1) = 1;
y2_explicit(1) = 1;
h = 0.1;
i = 1;
disp('Числа жёсткости') 
while(x_explicit(i) < 1) 
    g = [exp(x_explicit(i).^2), x_explicit(i); -1, 2];
    s = max(real(eig(g))) / min(real(eig(g)));
    fprintf('%d %d\n', i, s);
    i = i + 1;
    y1_explicit(i) = y1_explicit(i-1) + (y1_explicit(i-1) * exp(-x_explicit(i-1).^2) + x_explicit(i-1) * y2_explicit(i-1)) * h;
    y2_explicit(i) = y2_explicit(i-1) + (3 * x_explicit(i-1) - y1_explicit(i-1) + 2 * y2_explicit(i-1)) * h;
    x_explicit(i) = x_explicit(i-1) + h; 
end
fprintf('\n\n'); 
% Вывод найденных значений функции
fprintf('Значения функции y1 = f1(x) по явному методу Эйлера:\n\n');
for i = 1:length(x_explicit)
    disp([num2str(x_explicit(i)), '   ', num2str(y1_explicit(i))]);
    fprintf('\n');
end
fprintf('\n');
fprintf('Значения функции y2 = f2(x) по явному методу Эйлера:\n\n');
for i = 1:length(x_explicit)
    disp([num2str(x_explicit(i)), '   ', num2str(y2_explicit(i))]);
    fprintf('\n');
end
disp('___________________________________________________________________') 

% Эйлер неявный
x_implicit(1) = 0;
y1_implicit(1) = 1;
y2_implicit(1) = 1;
y1(1) = y1_implicit(1);
y2(1) = y2_implicit(1);
h = 0.1;
h2 = h / 2;
i = 1;
while(x_implicit(i) < 1)  
    i = i + 1;
    y1_implicit(i) = y1_implicit(i-1) + (y1_implicit(i-1) * exp(x_implicit(i-1).^2) + x_implicit(i-1) * y2_implicit(i-1)) * h;
    y2_implicit(i) = y2_implicit(i-1) + (3 * x_implicit(i-1) - y1_implicit(i-1) + 2 * y2_implicit(i-1)) * h;
    y1(i) = y1(i-1) + (y1(i-1) * exp(x_implicit(i-1).^2) + x_implicit(i-1) * y2(i-1) + y1(i-1) * exp(x_implicit(i-1).^2) + x_implicit(i-1) * y2(i-1)) * h2;
    y2(i) = y2(i-1) + (3 * x_implicit(i-1) - y1(i-1) + 2 * y2(i-1) + 3 * x_implicit(i-1) - y1_implicit(i-1) + 2 * y2_implicit(i-1)) * h2;
    x_implicit(i) = x_implicit(i-1) + h; 
end
% Вывод найденных значений функции 
fprintf('Значения функции y1 = f1(x) по неявному методу Эйлера:\n\n');
for i = 1:length(x_implicit)
    disp([num2str(x_implicit(i)), '   ', num2str(y1_implicit(i))]);
    fprintf('\n');
end
fprintf('Значения функции y2 = f2(x) по неявному методу Эйлера:\n\n');
for i = 1:length(x_implicit)
    disp([num2str(x_implicit(i)), '   ', num2str(y2_implicit(i))]);
    fprintf('\n');
end
disp('___________________________________________________________________')

% Стандартные методы ODE45
% Определение начальных условий
y0 = [1; 1];

% Решение системы дифференциальных уравнений
[x_ode45, y_ode45] = ode45(@myODEs, [0, 1], y0);

% Визуализация результатов
figure
hold on
plot(x_explicit, y1_explicit, 'r', 'DisplayName', 'Явный Эйлер y1');
plot(x_explicit, y2_explicit, 'g', 'DisplayName', 'Явный Эйлер y2');
plot(x_implicit, y1_implicit, 'b', 'DisplayName', 'Неявный Эйлер y1');
plot(x_implicit, y2_implicit, 'm', 'DisplayName', 'Неявный Эйлер y2');
plot(x_ode45, y_ode45(:,1), 'k--', 'DisplayName', 'MATLAB y1');
plot(x_ode45, y_ode45(:,2), 'c--', 'DisplayName', 'MATLAB y2');
grid on
xlabel('x');
ylabel('y');
legend;
title('Решение системы дифференциальных уравнений');
hold off

% Определение функции myODEs
function dydx = myODEs(x, y)
    dydx = [y(1) * exp(-x.^2) + x * y(2); 3 * x - y(1) + 2 * y(2)];
end