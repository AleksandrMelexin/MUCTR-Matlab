%��������� ������ ����������� �������� �������
%� ������� ������ �������� �������� �� �������
%���������������� �������
%������� ������� ������������
clear all
%���������� ������������ ������� ������� A=Q
nmax=200;
%������ �������� ������� A=Q
q=0.1;
%���������� ���� �������� ����������� ��������
%��� ������ ��������� �������� (�� 1 �� nmax)
for n=1:nmax
    %��������� ���������������� �������
    %A=Q ������� n
    A=zeros(n);
    for i=1:n
        A(i,i)=1+2*q;
    end
    for i=2:n
        A(i,i-1)=-q;
    end
    for i=1:(n-1)
        A(i,i+1)=-q;
    end
    %���� ������� ����������� �������� �������
    %�������� ��������
    for k=1:n
        %���������� ����������� ��������
        lambda(k)=1+4*q*sin((pi*k)/(2*(n+1)))^2;
        %��������� ����������� ��������
        lambda(k)=lambda(k)+1e-5*randn;
        %��������� ������ ����� ������� ���������
        %������ �������� ��������
        b=randn(n,1);
        %������������ ������������ ������� ������
        %�������� �������� � ����� ����������
        x0=b;
        for j=1:3
            %������ ������� �������� ��������
            x1=(A-lambda(k)*eye(n))\x0;
            x0=x1/norm(x1);
        end
        %����������� �������, � ������� ���������
        %�������� ����������� ������� ������� A=Q
        for j=1:n
            ev(j,k)=x0(j);
        end
    end
    %����������� �������� ��������� �����������
    %�������� ��� ������� ������� �������
    error(n)=0;
    for k=1:n
        error(n)=error(n)+...
                 norm((A-lambda(k)*eye(n))*ev(:,k));
    end
end
%�������� ����������� ����������� ����������
%����������� �������� �� ������� �������
plot(1:nmax,error);