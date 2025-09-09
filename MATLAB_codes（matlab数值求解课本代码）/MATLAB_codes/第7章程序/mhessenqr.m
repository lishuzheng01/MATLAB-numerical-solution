%上Hessenberg矩阵QR分解程序
function [Q,R]=mhessenqr(A)
%功能: 用Givens变换对上Hessenberg矩阵A进行QR分解
%输入: n阶上Hessenberg矩阵A, 其中A(i+1,i)=0, i>2.
%输出: 变换后的上Hessenberg形矩阵A.
n=size(A,1); Q=eye(n);
for i=1:n-1
    xi=A(i,i);  xk=A(i+1,i);
    if xk~=0
        d=sqrt(xi^2+xk^2);
        c=xi/d;  s=xk/d; J=[c, s;-s,c];
        A(i:i+1,i:n)=J*A(i:i+1,i:n);
        Q(1:n,i:i+1)=Q(1:n,i:i+1)*J';
    end
end
R=A;