%复化梯形公式MATLAB程序
function s=mtrap(f,a,b,n)
%用途: 用复化梯形公式求积分.
%格式: s=mtrap(f,a,b,n),  f为被积函数,a,b为积分区
%间的左右端点,n为区间的等分数,s返回积分近似值
h=(b-a)/n;
x=linspace(a,b,n+1);
y=feval(f,x);
s=0.5*h*(y(1)+2*sum(y(2:n))+y(n+1));