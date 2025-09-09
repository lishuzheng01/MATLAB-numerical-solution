%变步长梯形公式~MATLAB 程序
function [T,k,n]=mvtrap(f,a,b,epsi)
%用途:用变步长梯形公式求积分.
%格式:[T,k,n]=mvtrap(f,a,b,epsi),f为被积函数,a,b为积分区间的
%左右端点,epsi为控制精度,T返回积分近似值,k为对分区间的次数,
%n为区间的等分数
h=b-a;
T1=h*(feval(f,a)+feval(f,b))/2;
T=T1/2+(h/2)*feval(f,a+h/2);
k=1;n=2;
while(abs(T-T1)>epsi)
    h=h/2; T1=T;
    T=T1/2+(h/2)*sum(feval(f,a+h/2:h:b-h/2));
    k=k+1; n=2*n;
end