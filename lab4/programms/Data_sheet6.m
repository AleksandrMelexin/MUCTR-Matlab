%��������� ���������� ������� � ������� �����
%����������� ����� ������� ���������
%������� ������� ������������
clear all
%���������� ������������ ������� �������
nmax=403;
%���������� ���� �������������� ������ �������
%�������
for n=3:20:nmax
    %���������� ������ ������ ���������� ���
    %�������������� �������
    t0=cputime;
    %������� ��������� �������
    A=randn(n);
    %���������� ���� q-2 ��������
    %�������������� ���������
    for q=1:(n-2)
        u=0;
        for i=(q+1):n
            u=u+abs(A(i,q))^2;
        end
        %�������� ������ ��������� b(q+1,q)
        bq1q=sqrt(u);
        alpha=sqrt(2*bq1q*(bq1q+abs(A(q+1,q))));
        %��������� ������ ������� w � ���������
        %���������
        w=zeros(n,1);
        w(q+1)=(A(q+1,q)+bq1q*sign(A(q+1,q)))/alpha;
        for i=(q+2):n
            w(i)=A(i,q)/alpha;
        end
        %��������� ������� ��������� R ��� �������
        %�������� ������� w
        for i=1:n
            for j=1:n
                if i==j
                    R(i,j)=1-2*w(i)*w(j);
                else
                    R(i,j)=-2*w(i)*w(j);
                end
            end
        end
        %������������ �������� �������������� ������� A
        %� ������� ������� ��������� R
        B=R*A*R;
        A=B;
    end
    %���������� ����� ����������� ����������� �����������
    %�� �������������� ������� A � �������� �����
    %������������ ���� ������� ���������
    t(n)=cputime-t0;
    %����������� ��������� ������ ������� ������������
    %���������� �������� ��� ������������� ����������
    %����������� � ������� ������� hess
    t0=cputime;
    B=hess(randn(n));
    thess(n)=cputime-t0;
end
%������ � ����� ���� ������ ����������� ������� ������
%������������ ���������� �� �������������� ������� ��
%������� ������� ��� ����� ��������� (������ *) � ���
%���������� � MATLAB ��������� hess (����� "�����������")
semilogy(3:20:nmax,t([3:20:nmax]),'-*',...
         3:20:nmax,thess([3:20:nmax]),'-p');