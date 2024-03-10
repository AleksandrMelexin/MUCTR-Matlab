disp('задание 1')
A = [6, -1, -1; 1, -2, +3; +3, +4, +4]
b = [0, 1, -1]
A_inv = inv(A)
x = b*A_inv
n = cond(A)
err = x*A - b

disp('задание 2')
A = [9.1, 5.6, 7.8; 3.8, 5.1, 2.8; 4.1, 5.7, 1.2]
B = [9.8; 6.7; 5.8]
n = cond(A)
C = [A, B]
[a, b] = rref(C)
x = a(1:3, 4)
err = A*x - B

disp('задание 3')
A = [2.34, -1.42, -0.54, +0.21;
     1.44, -0.53, 1.43, -1.27;
     0.63, -1.32, -0.65, +1.43;
     0.54, 0.88, -0.67, -2.38]
f = [0.66; -1.44; 0.94; 0.73]
[L,U,P] = lu(A)
f = P*f
dim = rows(A)
%forward substitution
v = zeros(dim, 1)
for i = 1:dim
    v(i) = f(i)
    s = 0
    for j = 1:i-1
	s += L(i, j)*v(j)
    endfor
    v(i) -= s
    v(i) /= L(i, i)
endfor
%backward substitution
x = zeros(dim, 1)
for i = dim:-1:1
    x(i) = v(i)
    s = 0
    for j = dim:-1:i+1
	s += U(i, j)*x(j)
    endfor
    x(i) -= s
    x(i) /= U(i, i)
endfor
n = cond(A)
err = A*x - f

disp('задание 4')
A = [2, 1, 3, 0, 0, 0;
     0, 0, 3, 1, 1, 0;
     1, 0, 3, 0, 1, 0;
     0, 0, 1, 2, 0, 0;
     0, 1, 2, 0, 0, 0;
     0, 0, 1, 0, 0, 1;
     0, 0, 6, 0, 2, 1]
A_tr = transpose(A)
cond(A_tr(1:6, 2:6))
r = rank(A)
n = rows(A)
k = n - r
%алгоритм нахождения всех ненулевых подматриц
rows = 1:n
C = nchoosek(rows, k)
D = []
for i = 1:size(C)
    rows_to_delete = C(i, :)
    temp = A
    temp(rows_to_delete, :) = []
    [L, U, P] = lu(temp)
    values = diag(U)
    determinant = prod(values, 1)
    if abs(determinant) > 0.0001
	D = [D; rows_to_delete]
    endif
endfor
%загрузка пакета символик
echo off;
pkg load symbolic
syms Na2(CO3) HNO3 NaNO3 H2O CO2 CaO Ca(NO3_2)
M_all = [Na2(CO3), HNO3, NaNO3, H2O, CO2, CaO, Ca(NO3_2)]
equations = []
for l = 1:size(D)
    M_basis = [Na2(CO3), HNO3, NaNO3, H2O, CO2, CaO, Ca(NO3_2)]
    indexes = [D(l, 1); D(l, 2)]
    n = size(indexes)
    for i = 1:n
        M_basis(indexes(i)) = 0
    endfor
    for i = 1:n
        M_B = M_basis
        M_B(indexes(i)) = 1 
        res = zeros(1, columns(A))
        [dict] = solve(M_B * A == res)
        reagents = ""
        products = ""
        coef = 1
        for [val, key] = dict
            if (abs(val) < 1 && val > 0)
            disp(val)
                 coef = 1 / val
            endif
        endfor
        
        for [val, key] = dict
             val *= coef
            if (val == 0)
                continue
            endif
            if (val > 0)
                if(val == 1)
                    products = strcat(products," + ", char(key))
                else
                    products = strcat(products," + ", char(val),"*", char(key))
                endif
            else
                reagents = strcat(reagents," + ", char(-val),"*", char(key))
            endif
        endfor
        if(coef == 1)
            equation = strcat(reagents, " = ", products, " + ", char(M_all(indexes(i))))
        else
            equation = strcat(reagents, " = ", products, " + ", char(coef), "*", char(M_all(indexes(i))))
        endif
        equations = [equations; equation]
    endfor
endfor
