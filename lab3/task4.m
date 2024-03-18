% оперделяем матрицу соединений
% Na N C H O Ca
A = [2 0 1 0 3 0; 
    0 1 0 1 3 0; 
    1 1 0 0 3 0; 
    0 0 0 2 1 0;
    0 0 1 0 2 0;
    0 0 0 0 1 1; 
    0 2 0 0 6 1];
 
 
syms Na N C H O Ca
 
X = [Na; N; C; H; O; Ca];
 
disp(A*X);
mym = A*X;
 
syms CNa2O3 HNO3 NaO3 H2O CO2 CaO;
gh = [CNa2O3; HNO3; NaO3; H2O; CO2; CaO];

syms Na2CO3 HN03 NaNO3 H2O CaO CO2 CaN2O6
all_soed = [Na2CO3; HN03; NaNO3; H2O; CO2; CaO; CaN2O6];
%disp(gh); 

for i=1:6
    qw = mym(i,:);
    yt = gh(i);
    %fprintf('%s = %s', char(qw), char(yt));
    disp(' ');
end    
 
b = [  1 0 0 0 0;
       0 0 1 0 0;
       0 2 0 0 0;
       1 0 0 0 2;
       0 0 0 1 1;
    ];

b1 = [ 1 0 3 0 1 ;
       0 0 1 2 0 ;
       0 1 2 0 0 ;
       0 0 1 0 0 ;
       0 0 6 0 2 ;
    ];

M = [-2; -1; 0; 0; 0];

 
t = [ 1 0 0 0 0 0 0];
 
fprintf( '\nМатрица \n');
disp(b);
fprintf('________________________________________________________________\n');
fprintf( 'Определитель матрицы: %d\n',det(b)); 
fprintf( 'Норма матрицы: %d\n',norm(b)); 
fprintf( 'Число обусловленности матрицы: %d\n',cond(b)); 
fprintf( 'Вторая норма матриц: %d\n',norm(b,2)); 
fprintf( 'Число обусловленности матрицы для второй нормы %d\n',cond(b,2)); 
fprintf('________________________________________________________________\n\n');
 
x = inv(b)*M;
 
disp(t)
for i=1:5
    t(i+2) = x(i);
end
 
fprintf( 'Решение \n');
disp(t)

rez = t * all_soed;
%disp(rez);
fprintf('2*NaNC3 + CO2 + CaO = CaN2O6 + Na2CO3\n\n');

fprintf( 'Проверка условия B*A=0\n');
disp(t*A);
 
 
fprintf( '\n____________________________________________________________________________________________________________________\n\n');
 
b = [  1 0 0 0 0;
       0 0 1 0 0;
       0 2 0 0 0;
       1 0 0 0 2;
       0 0 0 1 1;
    ];
 
 
M = [0; 0; -1; -1; 0];
 
t = [ 0 1 0 0 0 0 0];

fprintf( '\nМатрица \n');
disp(b);
fprintf('________________________________________________________________\n');
fprintf( 'Определитель матрицы: %d\n',det(b)); 
fprintf( 'Норма матрицы: %d\n',norm(b)); 
fprintf( 'Число обусловленности матрицы: %d\n',cond(b)); 
fprintf( 'Вторая норма матриц: %d\n',norm(b,2)); 
fprintf( 'Число обусловленности матрицы для второй нормы %d\n',cond(b,2)); 
fprintf('________________________________________________________________\n\n');
 
 
 
 
x = inv(b) * M;
 
for i=1:5
    t(i+2) = x(i);
end
 
 
fprintf( 'Решение \n');
disp(t);

rez = t * all_soed;
%disp(rez);

reag = [H H];
prod = [H H];
fprintf('CaN2O6 + H2O2 = 2HNO3 + CaO\n\n');

fprintf( 'Проверка условия B*A=0\n');
disp(t*A);

