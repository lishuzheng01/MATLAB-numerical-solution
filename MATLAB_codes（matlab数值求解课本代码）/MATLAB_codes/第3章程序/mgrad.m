%梯度下降法MATLAB程序
function [x,iter]=mgrad(A,b,x,ep,N)
%用途：用梯度下降法解线性方程组Ax=b
%格式：[x,iter]=mgrad(A,b,x,ep,N)  其中A为系数矩阵, b为右端向
%量, x为初始向量(默认零向量), ep为精度(默认1e-6),N为最大迭代次
%数(默认1000次), 返回参数x,iter分别为近似解向量和迭代次数
if nargin<5, N=1000; end
if nargin<4, ep=1e-6; end
if nargin<3, x=zeros(size(b)); end
for iter=1:N
   r=b-A*x;
   if norm(r)<ep, break; end
   a=r'*r/(r'*A*r);
   x=x+a*r;
end