%幂法MATLAB程序
function [lam,v,k]=mepower(A,x,epsi,N)
%用途: 用幂法求矩阵的模最大特征值和对应的特征向量
%格式: [lam,v,k]=mepower(A,x,epsi,N)
%输入: A为n阶方阵, x为初始向量, epsi控制精度, N为最大迭代次数.
%输出: lam--按模最大的特征值, v--对应的特征向量, k--迭代次数.
if nargin<4, N=500; end
if nargin<3, epsi=1e-6; end
m=0; k=0; err=1;
while (k<N) & (err>epsi)
   v=A*x;
   [~,t]=max(abs(v));
   m1=v(t);  x=v/m1;
   err=abs(m1-m);
   m=m1;   k=k+1;
end
lam=m1; v=x;