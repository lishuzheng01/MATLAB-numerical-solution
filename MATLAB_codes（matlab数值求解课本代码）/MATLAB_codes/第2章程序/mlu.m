%LU分解MATLAB程序
function [x,L,U]=mlu(A,b)
%用途: 用LU分解法解方程组Ax=b
%格式: [x,l,u]=malu(A,b)  A为系数矩阵，b为右端向量
%返回: x--解向量，L--下三角矩阵，U--上三角矩阵
%LU分解
n=length(b); format short
U=zeros(n,n); L=eye(n,n);
U(1,:)=A(1,:);  L(2:n,1)=A(2:n,1)/U(1,1);
for k=2:n
    U(k,k:n)=A(k,k:n)-L(k,1:k-1)*U(1:k-1,k:n);
    L(k+1:n,k)=(A(k+1:n,k)-L(k+1:n,1:k-1)*U(1:k-1,k))/U(k,k);
end
%解下三角方程组Ly=b
y=zeros(n,1);
y(1)=b(1);
for k=2:n
    y(k)=b(k)-L(k,1:k-1)*y(1:k-1);
end
%解上三角方程组Ux=y
x=zeros(n,1);
x(n)=y(n)/U(n,n);
for k=n-1:-1:1
    x(k)=(y(k)-U(k,k+1:n)*x(k+1:n))/U(k,k);
end