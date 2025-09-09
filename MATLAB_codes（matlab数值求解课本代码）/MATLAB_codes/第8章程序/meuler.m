%改进的欧拉方法MATLAB程序
function [t,y]=meuler(df,tspan,y0,h)
%用途: 改进的欧拉方法解常微分方程y'=f(t,y), y(t0)=y0
%格式: [t,y]=meuler(df,tspan,y0,h), df为函数f(t,y), tspan为求解
%区间[t0,tn], y0为初值y(t0), h为步长, [t,y]返回节点和数值解矩阵
t=tspan(1):h:tspan(2);  y(1)=y0;
for n=1:(length(t)-1)
    k1=feval(df,t(n),y(n));
    y(n+1)=y(n)+h*k1;
    k2=feval(df,t(n+1),y(n+1));
    y(n+1)=y(n)+h*(k1+k2)/2;
end