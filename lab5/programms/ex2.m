function ex2
    % ������ ��������� f(x)=0,  ��� ��� f(x)= x^3 - cos(x) + 1 ������� ��������

    % ����� ������� f(x)
    f = inline('x.^3 - cos(x) + 1');
    root1 = bisec(f, -0.6, -0.4)
    root2 = bisec(f, -0.2, 0.2)
    
%   ����� ��������
function center = bisec(f, left, right)
    % ������������ ����� �������� � ��������� 2 eps    
    while right - left > eps * 2
        center = (right - left) / 2 + left;
        if f(center) * f(left) > 0
            left = center;
        else
            right = center;
        end
    end