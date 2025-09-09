%四阶经典龙格-库塔公式MATLAB程序
function [t,y]=m4rkutta(df,tspan,y0,h)
%用途: 4阶经典龙格-库塔公式解常微分方程y'=f(t,y), y(t0)=y0
%格式: [t,y]=m4rkutta(df,tspan,y0,h), df为函数f(t,y)表达式, 
%tspan为求解区间[t0,tn], y0为初值, h为步长, t为节点, y为数值解
t=tspan(1):h:tspan(2);  y(1)=y0;
for n=1:(length(t)-1)
    k1=feval(df,t(n),y(n));
    k2=feval(df,t(n)+h/2,y(n)+h/2*k1);
    k3=feval(df,t(n)+h/2,y(n)+h/2*k2);
    k4=feval(df,t(n+1),y(n)+h*k3);
    y(n+1)=y(n)+h*(k1+2*k2+2*k3+k4)/6;
end