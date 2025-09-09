%原点位移反幂法MATLAB程序
function [lam,v,k]=minvpower(A,x,alpha,epsi,N)
%用途: 用原点位移反幂法求矩阵的模最大特征值和对应的特征向量
%格式: [lam,v,k]=minvpower(A,x,alpha,eps,N)
%输入: A为n阶方阵,x为初始向量,epsi为控制精度(默认1.e-6),
%        N为最大迭代次数(默认500次),alpha为模最大的近似特征值.
%输出: lam--距离alpha最近的特征值,v--对应的特征向量,k--迭代次数
if nargin<5, N=500; end
if nargin<4,epsi=1e-6;end
m=0.5; k=0; err=1;
A=A-alpha*eye(length(x));
[L,U,P]=lu(A);
while (k<N) & (err>epsi)
   [~,t]=max(abs(x));
   m1=x(t); v=x/m1;
   z=L\(P*v);  x=U\z;
   err=abs(1/m1-1/m);
   k=k+1;  m=m1;
end
lam=alpha+1/m;