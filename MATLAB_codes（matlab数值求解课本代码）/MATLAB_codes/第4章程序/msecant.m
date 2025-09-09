%割线法MATLAB程序
function [x,iter]=msecant(f,x0,x1,ep,N)
%用途:用割线法求解非线性方程f(x)=0
%格式:[x,iter]=msecant(f,x0,x1,ep,N)  f为f(x)的表达式, x0, x1为
%迭代初值, ep为精度(默认1e-5), N为最大迭代次数(默认为500), 
%x,iter分别返回近似根和迭代次数
if nargin<5,N=500;end
if nargin<4,ep=1e-5;end
iter=0;
while iter<N
    x=x1-(x1-x0)*feval(f,x1)/(feval(f,x1)-feval(f,x0));
    if abs(x-x1)<ep,  break;  end
    x0=x1; x1=x; iter=iter+1;
end