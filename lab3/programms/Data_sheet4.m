%������� ����������� ������
n=30;
e=ones(n,1);
%������ ���������������� �������
A=spdiags([e e e],-1:1,n,n);
subplot(1,2,1); 
spy(A);
subplot(1,2,2);
%������ ��������� ������������
%����������� �������
spy(sprandsym(n,0.15));