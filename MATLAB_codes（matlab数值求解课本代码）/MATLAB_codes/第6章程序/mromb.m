%龙贝格积分公式MATLAB程序
function [T,n]=mromb(f,a,b,epsi)
%用途: 用龙贝格公式求积分
%格式: [T,n]=mromb(f,a,b,epsi), f是被积函数, [a,b]是积分
%区间,epsi控制精度, T返回积分近似值,n返回区间等分数
if nargin<4,epsi=1e-6;end
h=b-a;
R(1,1)=(h/2)*(feval(f,a)+feval(f,b));
n=1; J=0; err=1;
while (err>epsi)
    J=J+1; h=h/2; S=0;
    for i=1:n
        x=a+h*(2*i-1);
        S=S+feval(f,x);
    end
    R(J+1,1)=R(J,1)/2+h*S;
    for k=1:J
        R(J+1,k+1)=(4^k*R(J+1,k)-R(J,k))/(4^k-1);
    end
    err=abs(R(J+1,J+1)-R(J+1,J));
    n=2*n;
end
R;  %龙贝格积分表
T=R(J+1,J+1);