clc; clear;
syms x y(x)

% Задание дифференциального уравнения
eqn = x^3*(diff(y)-x) == y^2;

% Определение функций
y1 = (x^2/log(x))*(1-log(x));
y2 = x^2;

% Подстановка функций в уравнение и упрощение
eqn1 = subs(eqn, y(x), y1);
eqn2 = subs(eqn, y(x), y2);

% Проверка уравнений
sol1 = simplify(eqn1);
sol2 = simplify(eqn2);
if sol1 == symtrue
    disp('Первая функция является решением уравнения.')
else
    disp('Первая функция не является решением уравнения.')
end
disp(eqn1);
if sol2 == symtrue
    disp('Вторая функция является решением уравнения.')
else
    disp('Вторая функция не является решением уравнения.')
end
disp(eqn2);
