%矩阵上Hessenberg化的MATLAB程序
function A=mhessenb(A)
%用途: 用Householder变换化矩阵为上Hessenberg矩阵.
%输入: n阶实方阵A; 输出: A的上Hessenberg形
%调用函数: mhouseh.m
n = size(A,1);
for k=1:(n-2)
    x=A(k+1:n,k); H=mhouseh(x);
    A(k+1:n,1:n)=H*A(k+1:n,1:n);
    A(1:n,k+1:n)=A(1:n,k+1:n)*H;
end