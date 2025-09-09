%二分法MATLAB程序
function [x,iter]=mbisec(f,a,b,ep)
%用途: 用二分法求非线性方程f(x)=0有根区间[a,b]中的一个根
%格式: [x,iter]=mbisec(f,a,b,ep),  f为函数表达式，a,b为区间
%左右端点, ep为精度, x, iter分别返回近似根和二分次数
x=(a+b)/2.0;  iter=0;
while abs(feval(f,x))>ep|(b-a>ep)
    if feval(f,x)*feval(f,a)<0
        b=x;
    else
        a=x;
    end
    x=(a+b)/2.0; 
    iter=iter+1;
end