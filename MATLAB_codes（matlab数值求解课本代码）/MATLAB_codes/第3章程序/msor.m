%SOR迭代法MATLAB程序
function [x,iter]=msor(A,b,omega,x,ep,N)
%用途:用SOR迭代法解线性方程组Ax=b
%格式:[x,iter]=msor(A,b,omega,x,ep,N)  A为系数矩阵, b为右
%端向量, omega为松弛因子(默认1.05), x为初始向量(默认零向量),
%ep为精度(默认e-6), N为最大迭代次数(默认500次), 返回参数x,
%iter分别为近似解和迭代次数
if nargin<6, N=500; end
if nargin<5, ep=1e-6; end
if nargin<4, x=zeros(size(b)); end
if nargin<3, omega=1.05; end
D=diag(diag(A)); L=D-tril(A); U=D-triu(A);
for iter=1:N
   x=(D-omega*L)\(((1-omega)*D+omega*U)*x+omega*b);
   err=norm(b-A*x)/norm(b);
   if err<ep, break; end
end