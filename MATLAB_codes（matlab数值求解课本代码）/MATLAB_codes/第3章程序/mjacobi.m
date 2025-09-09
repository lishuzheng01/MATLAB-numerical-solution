%雅可比迭代法MATLAB程序
function [x,iter]=mjacobi(A,b,x,ep,N)
%用途：用雅可比迭代法解线性方程组Ax=b
%格式：[x,iter]=mjacobi(A,b,x,ep,N)  A为系数矩阵, b为右端向
%量, x为初始向量(默认零向量), ep为精度(默认1e-6), N为最大迭
%代次数(默认500次), 返回参数x,iter分别为近似解向量和迭代次数
if nargin<5, N=500; end
if nargin<4, ep=1e-6; end
if nargin<3, x=zeros(size(b)); end
D=diag(diag(A));
for iter=1:N
   x=D\((D-A)*x+b);
   err=norm(b-A*x)/norm(b);
   if err<ep, break; end
end