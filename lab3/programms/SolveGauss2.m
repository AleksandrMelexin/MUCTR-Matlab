function [ X ] = SolveGauss2(A,B)
%SLAE solving, Gauss method
%   Transform to Gauss view, A - matrix, B - free coefficients vector, X -
%   vector of solutions

%    ��� ����, ����� ������ ����� ������������� �����������, �������� ���
%    "%;" �� ";" ����� Ctrl-F

if (nargin < 2)
    disp ('Few arguments was specified, 2 matrices are required');
    %����� ���������� 2, ���� ������ - ����� ���������, ������� ������
    %�������������� �� �����
else
    [m,n] = size(A); % � m � n ���������� ������� ������� A - ����� ����� m � ����� �������� n
    if (m~=n)%���� ������� �� ����������
        if (m>n)
        disp('Equation system is overpredicted, please, reduce the number of  equations (strings)');
        %������� �������������, ����� �������� ������� ������ 0, ����� ��������� ����� ������, �.�. ���������
        else
        disp('Equation system has degrees of  freedom, please increase the number of equations (strings)'); 
        %������� ��������������, ����� �������� ������� ������ ����, �����
        %��������� ����� ������ (���������)
        end
        return; %����� �� �������
    end
    
    indexingfunc = @(Matrix,row,column) Matrix(row,column);      % inline function to index (take some element) from matrix
    %���������� � ��� �������  indexingfunc ��������� ��� ��������� -
    %������� (Matrix) � ��� ����� - ����� ������ (row) � ������� (column) � ���������� ������� Matrix(row,column)
    %��� ������� ��� ��� ������ �������� ������ ������� - ����� �������� -
    %����� ���������� ������� size(B), ��� ��� � Matlab ��� ��������
    %���������� � ���������� �������. ����� �������� (1,2) �� ������-������
    %size(B) ���� � ����� ����� �������� ������� B
    
     if (indexingfunc(size(B),1,2)>1) % size(B) returns two numbers
         %���� ����� �������� ������� B ������ 1
        disp('B must be vector (column)');
        %B ������ ���� ������-���������
        return;
     end
        
     A_00=[A B]  %;                                                                       %A_00 - ������� A � ���������� ������ �������� B
     
     %���������� ��� �����
     
     %���������� � ����������� �����
    X= zeros(m,1) %;                                                                     %������������� ������-������� ������� ������, ������ - m x 1
    for i = 1:m-1                                                                             %������ ���
        for k =i+1:m      
            %��� ������ k-�� ������, ������� � A ���� ��������������� i-��
            if (A(i,i)==0)
                [aa,bb]=exchange_rows(i,A,B)
                A=aa;
                B=bb;
            end
            elimination_coeff=A(k,i)./A(i,i) %;                               %��������� ����������� 
            A(k,1:n)=A(k,1:n)-A(i,1:n).*elimination_coeff   %;    %���������� ������� �� k-�� ������ i-��, ����������� �� elimination_coeff
           %k
           %i  %�������� ����� k � i
           B(k)=B(k)-B(i).*elimination_coeff  %;                         %������� �� k-��� �������� i-��, ����������� �� elimination_coeff
        end
    end
    
        %�������� � ����������� ����� ������� A_00
    for i = 1:m-1                                                                             
        for k =i+1:m                                                                         
                        %��� ������ k-�� ������, ������� � A ���� ��������������� i-��
            if (A_00(i,i)==0)
                [aa]=exchange_rows(i,A_00)
                A_00=aa;
            end
            elimination_coeff=A_00(k,i)./A_00(i,i) %;                                
            A_00(k,1:n+1)=A_00(k,1:n+1)-A_00(i,1:n+1).*elimination_coeff   %;    % �������� � ����������� �����
           %k
           %i  %�������� ����� k � i
        end
    end
     
    % ������� ���������-�������
    
    if (m_rank(A) == m_rank(A_00)) %���� ���� ������� A � ����������� ������� A ���������
        if  (m_rank(A)<n) %�� ������ ����� ����������� (����� ����� ��������)
            disp('Infinite count of solutions');  % ������� - ����������� �����
            return; %����� �� �������
        else    % ���� �� ����� ����� ����������� - ������� ������� ����� �� ������ ������ (�������� ���)
            for i=m:-1:1                                                                              %�������� ���
                X(i)=(B(i)-sum(A(i,i+1:m).*(X(i+1:m).')))./A(i,i) %;     % ������� � ��������� ������: A(m,m)*X(m)=B(m) => X(m) = B(m)/A(m,m)
            end                                                                                              %A(m-1,m-1)*X(m-1)+A(m-1,m)*X(m)=B(m-1) => X(m-1)=(B(m-1)-A(m-1,m)*X(m))/A(m-1,m-1)
        end
    else
            disp('No solutions'); 
            return; %����� �� �������
    end
        end                                                                                                  % ... A(i,i)*X(i)+A(i,i+1)*X(i+1)+...+A(i,m-1)*X(m-1)+A(i,m)*X(m)=B(i) =>
                                                                                                                %  X(i) = (B(i) - ����� �� p �� A(i,p)*X(p) )/A(i,i) - ����� �������
    end                                                                                                      %  ��� i ��������������� �������� ��� ������ ����� �����

function [ ret ] = m_rank(M) % ������� ������� ���� ������� M
[m,n] = size(M);
summ=0;
    for i=1:m
    summ=summ+(1-isequal(M(i,:),zeros(1,n))); %���� ������ ������� isequal ������ 1, 1-isequal=0
    end
    ret =summ %; %� ���������� summ ������������� ���� ������� M, ������� �� � ����� �� �������
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
