%复化中点公式MATLAB程序
function s=midpt(f,a,b,n)
%用途：用复化中点公式求积分.
%格式：s=mintm(f,a,b,n),f为被积函数,a,b为积分区
%间的左右端点,n为区间的等分数,s返回积分近似值
h=(b-a)/n;
x=linspace(a+h/2,b-h/2,n);
y=feval(f,x);
s=h*sum(y);