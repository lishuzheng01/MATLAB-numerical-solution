%基本QR方法MATLAB程序 
function [Iter,D]=mqrdec(A,epsi)
%用途: 用基本QR算法求实方阵的全部特征值.
%输入: n阶实方阵A, 控制精度epsi(默认是1.e-5)
%输出: 迭代次数Iter, A的全部特征值D
%调用函数: mhessenb.m, mhessenqr.m, eig--仅用于1,2矩阵
if nargin<2,epsi=1e-5;end
n=size(A,1); D=zeros(n,1); 
i=n; Iter=0;  %初始化
A=mhessenb(A);  %化矩阵A为Hessenberg形
%用基本QR算法进行迭代
while (1)
    if n<=2
        la=eig(A(1:n,1:n)); D(1:n)=la';  break;
    end
    Iter=Iter+1;
    [Q,R]=mhessenqr(A); %对上Hessenberg 矩阵做QR分解
    A=R*Q;  %做正交相似变换
    %下面的程序段是判断是否终止
    for k=n-1:-1:1
        if abs(A(k+1,k))<epsi
            if n-k<=2
                la=eig(A(k+1:n,k+1:n));
                j=i-n+k+1; D(j:i)=la';
                i=j-1; n=k; break;
            end
        end
    end
end