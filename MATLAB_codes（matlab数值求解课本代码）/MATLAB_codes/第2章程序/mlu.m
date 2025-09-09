%LU�ֽ�MATLAB����
function [x,L,U]=mlu(A,b)
%��;: ��LU�ֽⷨ�ⷽ����Ax=b
%��ʽ: [x,l,u]=malu(A,b)  AΪϵ������bΪ�Ҷ�����
%����: x--��������L--�����Ǿ���U--�����Ǿ���
%LU�ֽ�
n=length(b); format short
U=zeros(n,n); L=eye(n,n);
U(1,:)=A(1,:);  L(2:n,1)=A(2:n,1)/U(1,1);
for k=2:n
    U(k,k:n)=A(k,k:n)-L(k,1:k-1)*U(1:k-1,k:n);
    L(k+1:n,k)=(A(k+1:n,k)-L(k+1:n,1:k-1)*U(1:k-1,k))/U(k,k);
end
%�������Ƿ�����Ly=b
y=zeros(n,1);
y(1)=b(1);
for k=2:n
    y(k)=b(k)-L(k,1:k-1)*y(1:k-1);
end
%�������Ƿ�����Ux=y
x=zeros(n,1);
x(n)=y(n)/U(n,n);
for k=n-1:-1:1
    x(k)=(y(k)-U(k,k+1:n)*x(k+1:n))/U(k,k);
end