function ex3
    % Решить уравнение f(x)=0,  где где f(x)= x^3 - cos(x) + 1 методом Ньютона

    % Введём функцию f(x)
    f = inline('x.^3 - cos(x) + 1');
    % Её производная
    df = inline('3*x.^2 + sin(x)');
    root1 = newton(f, df, -0.5);
    % Проверим корни
    f(root1)
    root2 = newton(f, df, -0.1);
    % Проверим корни
    f(root2)    

%   Метод Ньютона
function root = newton(f, df, x0)
    root = x0 - f(x0) / df(x0);
    old_root = x0;
    while abs(old_root - root) > 2 * eps
        t = old_root;
        old_root = root;
        root = t - f(t) / df(t);
    end
    