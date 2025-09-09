%四阶亚当斯预报-校正”公式MATLAB 程序
function [t, y]=m4adams(df,tspan,y0,h)
%用途: 四阶亚当斯预报-校正公式解常微分方程y'=f(t, y), y(t0)=y0
%格式: [t, y]=m4adams(df,tspan,y0,h), df为函数f(t,y)表达式, 
%tspan为求解区间[t0,tn],y0为初值, h为步长, t为节点, y为数值解
t=tspan(1) : h : tspan(2);
[~,y]=m4rkutta(df,[t(1),t(4)],y0,h); %用四阶龙格-库塔公式计算4个初值
for n=4:(length(t)-1)
  p=y(n)+h/24*(55*feval(df,t(n),y(n))-59*feval(df,t(n-1),y(n-1)) ...
      +37*feval(df,t(n-2),y(n-2))-9*feval(df,t(n-3),y(n-3)));
  y(n+1)=y(n)+h/24*(feval(df,t(n-2),y(n-2))-5*feval(df,t(n-1),y(n-1))...
      +19*feval(df,t(n),y(n))+9*feval(df,t(n+1),p));
end