% задание 1
x = 2.5378;
Dx = 0.0001; % абсолютная погрешность x
dx = 3.94 * power(10, -5); % относительная погрешность x
y = 2.536;
Dy = 0.001; % абсолютная погрешность y
dy = 3.94 * power(10, -4); % относительная погрешность y
S1 = x + y; 
S2 = x - y;
% расчёт абсолютных погрешностей
DS1 = Dx + Dy;
DS2 = Dx - Dy;
% расчёт относительных погрешностей
dS1 = DS1/S1;
dS2 = DS2/abs(S2);
disp(['Предельная абсолютная погрешность суммы: ', num2str(DS1)]);
disp(['Предельная абсолютная погрешность разности: ', num2str(DS2)]);
disp(['Предельная относительная погрешность суммы: ', num2str(dS1)]);
disp(['Предельная относительная погрешность разности: ', num2str(dS2)]);