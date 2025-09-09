%复化辛普森公式MATLAB程序
function s=msimp(f,a,b,n)
%用途: 用复化辛普森公式求积分.
%格式: s=msimp(f,a,b,n)  f为被积函数;a,b为积分区
%间的左右端点,n为区间的等分数,s返回积分近似值
h=(b-a)/n;
x=linspace(a,b,2*n+1);
y=feval(f,x);
s1=sum(y(3:2:2*n-1));
s2=sum(y(2:2:2*n));
s=(h/6)*(y(1)+2*s1+4*s2+y(n+1));