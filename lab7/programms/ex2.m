% �������-�������� ������������ ������� �����
% ����� ������� �����
f = inline('1./(1+25*x.^2)');

% �������� ������� ��������
x = linspace(-1, 1, 10);
y = f(x);

% �������� ������ �������-�������� �������������
plot(x, y);
