function ex3
    % ���������� ��������� ����� f(x),  ��� ��� f(x)= (tan(x)^2 -1)^2 ������� �������

    % ����� ������� f(x)
    syms x;
    f = (tan(x)^2 -1)^2;
    [root1, iter1] = newton(f, 0.5, 1);
    [root2, iter2] = newton(f, 0.5, 2);
    [root3, iter3] = newton(f, 0.5, 3);
    % �� m, ��� ������ �������� - ��������� �����
    iter1
    iter2
    iter3
    % ������ ����� �������� ��� m = 2, ������ ��������� ����� - 2

%   ����� �������
function [root, iter] = newton(f, x0, m)
    df = diff(f);
    root = x0 - subs(f, x0) / subs(df, x0);
    old_root = x0;
    iter = 0;
    while abs(subs(f, old_root)) > 10 * eps
        t = old_root;
        old_root = root;
        root = t - m*subs(f, t) / subs(df, t);
        iter = iter + 1;
        if iter > 200
            warning('Max number of iterations reached');
            root = NaN;
            return;
        end
    end
