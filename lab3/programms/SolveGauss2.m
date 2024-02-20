function [ X ] = SolveGauss2(A,B)
%SLAE solving, Gauss method
%   Transform to Gauss view, A - matrix, B - free coefficients vector, X -
%   vector of solutions

%    Для того, чтобы убрать вывод промежуточных результатов, замените все
%    "%;" на ";" через Ctrl-F

if (nargin < 2)
    disp ('Few arguments was specified, 2 matrices are required');
    %число аргументов 2, если меньше - вывод сообщения, функция дальше
    %рассчитываться не будет
else
    [m,n] = size(A); % в m и n помещаются размеры матрицы A - число строк m и число столбцов n
    if (m~=n)%если матрица не квадратная
        if (m>n)
        disp('Equation system is overpredicted, please, reduce the number of  equations (strings)');
        %Система переопрделена, число степеней свободы меньше 0, нужно уменьшить число связей, т.е. уравнений
        else
        disp('Equation system has degrees of  freedom, please increase the number of equations (strings)'); 
        %Система недоопределена, число степеней свободы больше нуля, нужно
        %увеличить число связей (уравнений)
        end
        return; %выход из функции
    end
    
    indexingfunc = @(Matrix,row,column) Matrix(row,column);      % inline function to index (take some element) from matrix
    %встроенная в код функция  indexingfunc принимает три аргумента -
    %матрицу (Matrix) и два числа - номер строки (row) и столбца (column) и возвращает элемент Matrix(row,column)
    %это сделано так как нельзя получить второй элемент - число столбцов -
    %после выполнения функции size(B), так как в Matlab нет операции
    %индексации к результату функции. Выбор элемента (1,2) из вектор-строки
    %size(B) даст в итоге число столбцов матрицы B
    
     if (indexingfunc(size(B),1,2)>1) % size(B) returns two numbers
         %если число столбцов матрицы B больше 1
        disp('B must be vector (column)');
        %B должен быть вектор-стоолбцом
        return;
     end
        
     A_00=[A B]  %;                                                                       %A_00 - матрица A с дописанным справа столбцом B
     
     %собственно сам метод
     
     %приведение к ступенчатой форме
    X= zeros(m,1) %;                                                                     %инициализация вектор-столбца решений нулями, размер - m x 1
    for i = 1:m-1                                                                             %прямой ход
        for k =i+1:m      
            %для каждой k-ой строки, лежащей в A ниже рассматриваемой i-ой
            if (A(i,i)==0)
                [aa,bb]=exchange_rows(i,A,B)
                A=aa;
                B=bb;
            end
            elimination_coeff=A(k,i)./A(i,i) %;                               %вычислить коэффициент 
            A(k,1:n)=A(k,1:n)-A(i,1:n).*elimination_coeff   %;    %поэлементо вычесть из k-ой строки i-ую, домноженную на elimination_coeff
           %k
           %i  %тестовый вывод k и i
           B(k)=B(k)-B(i).*elimination_coeff  %;                         %вычесть из k-ого элемента i-ый, домноженный на elimination_coeff
        end
    end
    
        %приведем к ступенчатой форме матрицу A_00
    for i = 1:m-1                                                                             
        for k =i+1:m                                                                         
                        %для каждой k-ой строки, лежащей в A ниже рассматриваемой i-ой
            if (A_00(i,i)==0)
                [aa]=exchange_rows(i,A_00)
                A_00=aa;
            end
            elimination_coeff=A_00(k,i)./A_00(i,i) %;                                
            A_00(k,1:n+1)=A_00(k,1:n+1)-A_00(i,1:n+1).*elimination_coeff   %;    % приводим к ступенчатой форме
           %k
           %i  %тестовый вывод k и i
        end
    end
     
    % теорема Кронекера-Капелли
    
    if (m_rank(A) == m_rank(A_00)) %если ранг матрицы A и расширенной матрицы A совпадают
        if  (m_rank(A)<n) %но меньше числа неизвестных (равно числу столбцов)
            disp('Infinite count of solutions');  % решений - бексонечное число
            return; %выход из функции
        else    % если же равно числу неизвестных - решения находим далее по методу Гаусса (обратный ход)
            for i=m:-1:1                                                                              %обратный ход
                X(i)=(B(i)-sum(A(i,i+1:m).*(X(i+1:m).')))./A(i,i) %;     % начиная с последней строки: A(m,m)*X(m)=B(m) => X(m) = B(m)/A(m,m)
            end                                                                                              %A(m-1,m-1)*X(m-1)+A(m-1,m)*X(m)=B(m-1) => X(m-1)=(B(m-1)-A(m-1,m)*X(m))/A(m-1,m-1)
        end
    else
            disp('No solutions'); 
            return; %выход из функции
    end
        end                                                                                                  % ... A(i,i)*X(i)+A(i,i+1)*X(i+1)+...+A(i,m-1)*X(m-1)+A(i,m)*X(m)=B(i) =>
                                                                                                                %  X(i) = (B(i) - Сумма по p от A(i,p)*X(p) )/A(i,i) - общая формула
    end                                                                                                      %  где i последовательно проходит все строки снизу вверх

function [ ret ] = m_rank(M) % функция находит ранг матрицы M
[m,n] = size(M);
summ=0;
    for i=1:m
    summ=summ+(1-isequal(M(i,:),zeros(1,n))); %если строка нулевая isequal выдаст 1, 1-isequal=0
    end
    ret =summ %; %в переменной summ сформировался ранг матрицы M, который мы и вернём из функции
end

function [X,varargout]=exchange_rows(i,A,varargin)
l=size(A,1);
B=[];
    for j=i:l
        if (A(j,i)~=0)
            kk=A(j,:);
            A(j,:)=A(i,:);
            A(i,:)=kk;
            if (nargin==3)
                B=varargin{1};
                kb=B(j);
                B(j)=B(i);
                B(i)=kb;
            end
            break;
        end

    end
    X=A;
    if (nargin==3)
        varargout{1}=B;
    end
end
