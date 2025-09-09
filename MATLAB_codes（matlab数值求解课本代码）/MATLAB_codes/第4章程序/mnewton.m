%牛顿法MATLAB程序
function [x,iter]=mnewton(f,df,x0,ep,N)
%用途: 用牛顿法求解非线性方程f(x)=0
%格式: [x,iter]=mnewton(f,df,x0,ep,N)
%f和df分别为表示f(x)及其导数, x0为迭代初值, ep为精度(默认1e-5), 
%N为最大迭代次数(默认为500), x,iter分别返回近似根和迭代次数
if nargin<5,N=500;end
if nargin<4,ep=1e-5;end
iter=0;
while iter<N
    x=x0-feval(f,x0)/feval(df,x0);
    if abs(x-x0)<ep
        break;
    end
    x0=x; iter=iter+1;
end