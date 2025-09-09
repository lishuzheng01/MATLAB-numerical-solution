%列主元LU分解法MATLAB程序
function [x,L,U,P]=mzlu(A,b)
%用途: 用列主元LU分解法解方程组Ax=b
%格式: [x,L,U,P]=mzlu(A,b),  A为系数矩阵, b为右端向量,
%返回: x-解向量, L-单位下三角矩阵, U-上三角矩阵,
%P-选主元时记录行交换的置换阵
n=length(b); 
P=eye(n); %P记录选择主元时候所进行的行变换
%列主元LU分解
for k=1:n
    A(k:n,k)=A(k:n,k)-A(k:n,1:k-1)*A(1:k-1,k);
    [s,m]=max(abs(A(k:n,k)));   %选列主元
    m=m+k-1;
    if m~=k
        A([k m],:)=A([m k],:);
        P([k m],:)=P([m k],:);
    end
    A(k+1:n,k)=A(k+1:n,k)/A(k,k);
    A(k,k+1:n)=A(k,k+1:n)-A(k,1:k-1)*A(1:k-1,k+1:n);
end
L=tril(A,-1)+eye(n,n); U=triu(A);
%解单位下三角矩阵 Ly=b
newb=P*b;  y=zeros(n,1);
for k=1: n
    j=1: k-1;
    y(k)=newb(k)-L(k,j)*y(j);
end 
%解上三角方程组Ux=y
x=zeros(n,1);
for k=n:-1:1
     j=k+1:n;
     x(k)=(y(k)-U(k,j)*x(j))/U(k,k);
end