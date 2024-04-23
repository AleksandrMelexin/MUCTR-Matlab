clear; clc;
mass_heatCapacity = [2.799505654 2.865976811 2.936424531 3.010995579 3.089799019 3.172890855 3.260255665 3.351785366 3.44725572 3.546301795 3.648394432 3.752820666 3.858671918 3.964844364 4.070055767 4.172881981 4.27181405 4.365333422 4.451998938 4.530535861 4.599915455];
pressure = [4865 5065 5265 5465 5665 5865 6065 6265 6465 6665 6865 7065 7265 7465 7665 7865 8065 8265 8465 8665 8865];
weights = [0.5 0.2 0.9 0.5 0.8 0.8 0.2 0.5 0.4 0.2 0.6 0.5 0.2 0.7 0.1 0.8 0.6 0.9 0.6 0.9 0.4];
% 1 
% ���������� ������� �������� ���������
dt = diff(mass_heatCapacity);
disp(['������� �������� ���������: ',dt]);
disp(dt);
diff_table = diffTable(pressure, mass_heatCapacity, transpose(weights));  % ����� ������� diffTable ��� ���������� ������� �������� ���������

% ������������ ������� �������� (n-2), ��� n - ���������� ������� �����
max_degree = length(pressure) - 2;

% ����������� ������� �������� �� ������ ������������ �������� � �������
degree = find(max(abs(diff_table)) > 1e-6, 1, 'last');
degree = min(degree, max_degree);  % ������������ ������� �������� ������������ ���������

% ����� ������� ��������
disp(['������������� ������� ��������: ', num2str(degree)]);

% 2 
% ���������� ������������������ ��������
% ��� ����� ������� ������������� � �������������� 
% ������������ �����������

% �������� ������� �����������
W = vander(pressure);
W = W(:,  end-4:end); % ����� ��������� ��������, ��������������� �������� ��������

% ������� ������� �������� ���������
coefficients = W \ mass_heatCapacity';

% ���������� �������� �������� ��� ���������� �������
p_density_w = polyval(flip(coefficients'), pressure);

% ���������� �������
figure;
plot(pressure, mass_heatCapacity, 'o', pressure, p_density_w, '-');
xlabel('�������� � ���');
ylabel('�����������');
title(['����������������� �������(����������) ������� ', num2str(degree)]);
legend('����������������� ������', '����������������� ������� (����������)');
grid on;

% 3 
% ��������� ����������������� ������� 
% ��� ����� ������� ������������� � �������������� 
% ����������� ���������� MATLAB

% ������������� ���������
p = polyfit(pressure, mass_heatCapacity, degree);
% ���������� �������� �������� ��� ���������� �������
p_density_s = polyval(p, pressure);
% ���������� �������
figure;
plot(pressure, mass_heatCapacity, 'o', pressure, p_density_s, 'r-');
xlabel('�������� � ���');
ylabel('�����������');
title(['����������������� ������� ������� ', num2str(degree)]);
legend('����������������� ������', '����������������� ������� ��������');
grid on;

% 4 
% ��������� ����������������� ������� 
% � ������ ������� ������������� � �������������� 
% ������� spap2

% �������� ������� �������
sp = spap2(pressure, 1, pressure, mass_heatCapacity, weights);

% ������ �������� ��������� ��� ���������� �������
p_density2 = fnval(sp, pressure);

% ���������� �������
figure;
plot(pressure, mass_heatCapacity, 'o', pressure, p_density2, 'g-');
xlabel('�������� � ���');
ylabel('�����������');
title(['����������������� �������(spap2) ������� ', num2str(degree)]);
legend('����������������� ������', '����������������� ������� spap2');
grid on;

% 5 ��������� ����������������� ������� 
% � ������ ������� ������������� � �������������� 
% ������� fminsearch

% ������� ��� ����������� (����� ��������� ���������)
fun = @(coefficients) sum(weights .* (polyval(fliplr(coefficients'), pressure) - mass_heatCapacity).^2);

% ��������� ����������� ��� ������������� ��������
initial_guess = rand(1, degree+1); % ��������� ��������� �����������

% ����������� ������� � �������������� fminsearch
optimal_coefficients = fminsearch(fun, initial_guess);

% ���������� �������� �������� ��� ���������� �������
p_density_f = polyval(flip(optimal_coefficients'), pressure);

% ���������� �������
figure;
plot(pressure, mass_heatCapacity, 'o', pressure, p_density_f, 'b-');
xlabel('�������� � ���');
ylabel('�����������');
title(['����������������� �������(fminsearch) ������� ', num2str(degree)]);
legend('����������������� ������', '����������������� ������� fminsearch');
grid on;

% 6 
% ������� �������� �������������
% ���������� MSE ��� ������� ������
mse_polyfit = mean((polyval(p, pressure) - mass_heatCapacity).^2);
mse_vander = mean((polyval(flip(coefficients'), pressure) - mass_heatCapacity).^2);
mse_spap2 = mean((fnval(sp, pressure) - mass_heatCapacity).^2);
mse_fminsearch = mean((polyval(flip(optimal_coefficients'), pressure) - mass_heatCapacity).^2);

% ����� �����������
disp('������������������ ������ (MSE) ��� ������� ������:');
disp(['1. ������� ������� ���������� ��������� (��� ����� ������� �������������): ', num2str(mse_polyfit)]);
disp(['2. ������� � �������������� ������������ �����������: ', num2str(mse_vander)]);
disp(['3. ������� � �������������� ������� spap2: ', num2str(mse_spap2)]);
disp(['4. ������������� � �������������� ������� fminsearch: ', num2str(mse_fminsearch)]);

% 1 � 2, 14 � 15 �� 3 ������
% ���������� ���������� ����������� ��� ������� 1 � 2
abs_error_1_2_polyfit = abs(polyval(p, pressure(2)) - mass_heatCapacity(2));
abs_error_1_2_vander = abs(polyval(flip(coefficients'), pressure(2)) - mass_heatCapacity(2));

% ���������� ���������� ����������� ��� ������� 14 � 15
abs_error_14_15_polyfit = abs(polyval(p, pressure(15)) - mass_heatCapacity(15));
abs_error_14_15_vander = abs(polyval(flip(coefficients'), pressure(15)) - mass_heatCapacity(15));

% ���������� �������� ����������� �������� ����������� ��� ������� ������
mean_abs_error_1_2_polyfit = mean(abs_error_1_2_polyfit);
mean_abs_error_1_2_vander = mean(abs_error_1_2_vander);
mean_abs_error_14_15_polyfit = mean(abs_error_14_15_polyfit);
mean_abs_error_14_15_vander = mean(abs_error_14_15_vander);

% ����� �����������
disp('������� ���������� �������� ����������� ��� ������� ������:');
disp(['1. ������� ������� ���������� ��������� (��� ����� ������� �������������), ����� 1 � 2 �������: ', num2str(mean_abs_error_1_2_polyfit)]);
disp(['2. ������� � �������������� ������������ �����������, ����� 1 � 2 �������: ', num2str(mean_abs_error_1_2_vander)]);
disp(['3. ������� ������� ���������� ��������� (��� ����� ������� �������������), ����� 14 � 15 �������: ', num2str(mean_abs_error_14_15_polyfit)]);
disp(['4. ������� � �������������� ������������ �����������, ����� 14 � 15 �������: ', num2str(mean_abs_error_14_15_vander)]);

% 7 
% ��������� �� ����� ������� ���������� ������� 
% � ����������� �� ��� �������� ������� � ���� ��������, ������� �������, 
% ������� �������.

figure;
hold on
grid on
% ����� � ����������
plot(pressure, mass_heatCapacity, 'o', pressure, p_density_w, '-');
% ��������
plot(pressure, p_density_s, '-');
% ����2
plot(pressure, p_density2, '-');
% ��������
plot(pressure, p_density_f, '-');
xlabel('�������� � ���');
ylabel('�����������');
title('����������������� ��������');
legend('����������������� ������', '������� + ����������', '����������� ������', 'spap2', 'fminsearch');
hold off

% 8 
% ��������� �������, ���������������� ����������������� ������, 
% �� � ���� ��������, � ���� ������ ������� 
% � ������ ������� ������������� � �������������� ������� fminsearch

% ������� � �����������(���)
obj_function = @(params) sum(weights .* (mass_heatCapacity - (params(1) .* pressure + params(2)) ./ (params(3) .* pressure + params(4))).^2);

% �������
initial_guess = [1, 1, 1, 1];

% ����������� fminsearch
params = fminsearch(obj_function, initial_guess);

% ������������� ����������
a = params(1);
b = params(2);
c = params(3);
d = params(4);

% ����� ����������� ����������
disp(['Optimized Parameters: a = ', num2str(a), ', b = ', num2str(b), ', c = ', num2str(c), ', d = ', num2str(d)]);

% ������ ������ � ������
p_range = linspace(min(pressure), max(pressure), 100);
fitted_density = (a .* p_range + b) ./ (c .* p_range + d);

figure;
plot(pressure, mass_heatCapacity, 'bo', p_range, fitted_density, 'm-');
grid on
xlabel('�������� (���)');
ylabel('�����������');
legend('����������������� ������', '������');
title('����������������� ������� ��� ��������');

% 9 
% ��������� �������, ���������������� ����������������� ������, 
% � ���� �������� ��������, � ������ ������� �������������

% ��������� ������� �������� ��������
degree = 5;

% ��������� ���� 
x = 2 * ((pressure - min(pressure)) / (max(pressure) - min(pressure))) - 1;

% ���������� ���������� � ������� ��������������� ���
p_coefficients = polyfit(x, mass_heatCapacity, degree);

% ��������� � �������� ������ ������
fitted_density = polyval(p_coefficients, x);

% ���������� �������� ������ � ����������� ������
figure;
plot(pressure, mass_heatCapacity, 'bo', pressure, fitted_density, 'c-');
grid on
xlabel('�������� (���)');
ylabel('�����������');
legend('����������������� ������', '�������');
title('����������������� ������� ��������');


function diff_table = diffTable(x, y, w)
    % ������������� ������� �������� ���������
    n = length(x);
    diff_table = zeros(n, n);
    diff_table(:, 1) = w .* y';

    % ���������� �������� ���������
    for j = 2:n
        for i = 1:n-j+1
            w_factor = w(i) / w(i+j-1);  % ����������� ����� ��� ������� ��������
            diff_table(i, j) = w_factor * (diff_table(i+1, j-1) - diff_table(i, j-1));
        end
    end
end
