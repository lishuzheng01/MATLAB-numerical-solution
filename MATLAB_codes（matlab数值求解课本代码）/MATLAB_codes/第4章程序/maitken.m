%艾特金-史蒂文森加速法MATLAB程序
function [x,iter]=maitken(phi,x0,ep,N)
%用途: 用艾特金-史蒂文森加速法求f(x)=0的解
%格式: [x,iter]=maitken(phi,x0,ep,N), phi为迭代函数, x0为迭代初值,ep为精
%度(默认1e-5), N为最大迭代次数(默认500), x,iter分别返回近似根和迭代次数
if nargin<4,N=500;end
if nargin<3,ep=1e-5;end
iter=0;
while iter<N
   y=feval(phi,x0);
   z=feval(phi,y);
   x=x0-(y-x0)^2/(z-2*y+x0);
   if abs(x-x0)<ep, break;  end
   x0=x; iter=iter+1;
end