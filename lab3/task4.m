% оперделяем матрицу соединений
% Na N C H O Ca
A = [2 0 1 0 3 0; 
    0 1 0 1 3 0; 
    1 1 0 0 3 0; 
    0 0 0 2 1 0;
    0 0 1 0 2 0;
    0 0 0 0 1 1; 
    0 2 0 0 6 1];
A_tr = transpose(A)
cond(A_tr(1:6, 2:6))
r = rank(A)
n = 7;
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
    end
end
%загрузка пакета символик
echo off;
pkg load symbolic;
syms Na2(CO3) HNO3 NaNO3 H2O CO2 CaO Ca(NO3_2)
M_all = [Na2(CO3), HNO3, NaNO3, H2O, CO2, CaO, Ca(NO3_2)]
equations = []
for l = 1:size(D)
    M_basis = [Na2(CO3), HNO3, NaNO3, H2O, CO2, CaO, Ca(NO3_2)]
    indexes = [D(l, 1); D(l, 2)]
    n = size(indexes)
    for i = 1:n
        M_basis(indexes(i)) = 0
    end
    for i = 1:n
        M_B = M_basis
        M_B(indexes(i)) = 1 
        res = zeros(1, columns(A))
        [dict] = solve(M_B * A == res)
        reagents = ""
        products = ""
        coef = 1
        for val = key:dict
            if (abs(val) < 1 && val > 0)
            disp(val)
                 coef = 1 / val
            end
        end
        
        for val = key:dict
             val *= coef
            if (val == 0)
                continue
            end
            if (val > 0)
                if(val == 1)
                    products = strcat(products," + ", char(key))
                else
                    products = strcat(products," + ", char(val),"*", char(key))
                end
            else
                reagents = strcat(reagents," + ", char(-val),"*", char(key))
            end
        end
        if(coef == 1)
            equation = strcat(reagents, " = ", products, " + ", char(M_all(indexes(i))))
        else
            equation = strcat(reagents, " = ", products, " + ", char(coef), "*", char(M_all(indexes(i))))
        end
        equations = [equations; equation]
    end
end
