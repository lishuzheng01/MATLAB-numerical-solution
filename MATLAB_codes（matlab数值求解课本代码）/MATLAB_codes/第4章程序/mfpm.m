%不动点迭代法MATLAB程序
function [x,iter]=mfpm(phi,x0,ep,N)
%用途: 用不动点迭代法求非线性方程f(x)=0有根区间[a,b]中的一个根
%格式: [x,iter]=mfpm(phi,x0,ep,N), phi为迭代函数, x0为初值, ep为精度
%(默认1e-5), N为最大迭代次数(默认500), x,iter分别为近似根和迭代次数
if nargin<4 N=500;end
if nargin<3 ep=1e-6;end
iter=0;
while iter<N
   x=feval(phi,x0);
   if abs(x-x0)<ep
      break;
   end
   x0=x; iter=iter+1;
end