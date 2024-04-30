clc; clear;
f = @(x) x.^2 + exp(x+3); % исхожная функция
f2 = @(x) 2 + exp(x+3); % 2-ая производная
f4 = @(x) exp(x+3); % 4-ая производная
% Интервал интегрирования
a = -5;
b = 5;

tol = 0.01; % Точность вычисления интеграла

% поиск максимального значения 2-ой производной
m2 = NaN;
for i=a:0.1:b
    if isnan(m2)
        m2 = f2(i);
    end
    if f2(i) > m2
        m2 = f2(i);
    end    
end

h = sqrt((12 * tol)/(b-a)*m2); % расчёт шага h по формуле 5
fprintf('Шаг h для метода трапеций: %.10f\n', h);
fprintf('Максимальное значение 2-ой производной на интервале (a; b): %.10f\n', m2);
I_trapezoidal = (h/2) * (f(a) + f(b)); % Метод трапеций с уточнением формулы по процедуре Рунге
n = 1;
err = tol + 1;
while err > tol
    h = h/2;
    x = a + h : h : b - h;
    I_new = (h/2) * (f(a) + 2 * sum(f(x)) + f(b));
    err = abs(I_new - I_trapezoidal);
    I_trapezoidal = I_new;
    n = n + 1;
end


% График зависимости точности от шага h
h_values = zeros(1, n);
err_values = zeros(1, n);
h = (b - a);
for i = 1 : n
    h_values(i) = h;
    x = a + h : h : b - h;
    I_new = (h/2) * (f(a) + 2 * sum(f(x)) + f(b));
    err_values(i) = abs(I_new - I_trapezoidal);
    h = h/2;
end

figure;
plot(h_values, err_values);
xlabel('Шаг h');
ylabel('Погрешность');
title('Зависимость точности интегрирования методом трапеций от шага h');

% поиск максимального значения 2-ой производной
m4 = NaN;
for i=a:0.1:b
    if isnan(m4)
        m4 = f4(i);
    end
    if f4(i) > m4
        m4 = f4(i);
    end    
end

h = sqrt((180 * tol)/(b-a)*m4); % расчёт шага h по формуле 5
fprintf('Шаг h для метода Симпсона: %.10f\n', h);
fprintf('Максимальное значение 4-ой производной на интервале (a; b): %.10f\n', m4);
I_simpson = h/3 * (f(a) + 4*f((a+b)/2) + f(b)); % метод Симпсона
n = 1;
err = tol + 1;
while err > tol
    h = h/2;
    x = a + h : 2* h : b - h;
    I_new = h/3 * (f(a) + 4 * sum(f(x)) + 2 * sum(f(x + h)) + f(b));
    err = abs(I_new - I_simpson);
    I_simpson = I_new;
    n = n + 1;
end

I_matlab = integral(f, a, b); % Использование стандартных функций MATLAB для вычисления интеграла

% Вывод результатов
fprintf('Интеграл методом трапеций: %.10f\n', I_trapezoidal);
fprintf('Интеграл методом Симпсона: %.10f\n', I_simpson);
fprintf('Интеграл, вычисленный с помощью стандартных функций MATLAB: %.10f\n', I_matlab);

%вычисление неопределённого и несобственного интегралов
syms x a p;
f_2 = (a^x) * exp(-x);
f_3 = (1+x)/((x+a).^(p+1));
res = int(f_2, x);
fprintf('Неопредеённый интеграл: %s\n', res);
res = int(f_3, x, 0, Inf);
fprintf('Несобственный интеграл: %s\n', res);
