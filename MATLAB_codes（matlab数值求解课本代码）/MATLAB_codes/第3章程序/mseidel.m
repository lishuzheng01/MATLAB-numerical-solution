%高斯-赛德尔迭代法MATLAB程序
function [x,iter]=mseidel(A,b,x,ep,N)
%用途：用高斯-赛德尔迭代法解线性方程组Ax=b
%格式：[x,iter]=mseidel(A,b,x,ep,N)  A为系数矩阵, b为右端向
%量, x为初始向量(默认零向量), ep为精度(默认1e-6), N为最大迭
%代次数(默认500次), 返回参数x,iter分别为近似解向量和迭代次数
if nargin<5, N=500; end
if nargin<4, ep=1e-6; end
if nargin<3, x=zeros(size(b)); end
D=diag(diag(A)); L=D-tril(A); U=D-triu(A);
for iter=1:N
   x=(D-L)\(U*x+b);
   err=norm(b-A*x)/norm(b);
   if err<ep, break; end
end