function U=finitedif(f,g,c,a,b,n,m)
%功能: 用差分法求解波动方程u_tt=c^2*u_xx
%输入: f=u(x,0), 初始条件函数,f为函数句柄;
%        g=u_t(x,t),导数初始条件,g为函数句柄;
%        a,b分别为空间域[0,a]和时间域[0,b]右端点
%        n,m分别为[0,a]和[0,b]上的等分节点数
%输出: U为解矩阵,行为空间节点,列为时间节点
%初始化
h=a/(n-1); k=b/(m-1);
r=c^2*k^2/h^2;  U=zeros(n,m);
%计算第一层和第二层
for i=2:n-1
    U(i,1)=feval(f,h*(i-1)); %第一层
    U(i,2)=(1-r)*feval(f,h*(i-1))+k*feval(g,h*(i-1)) ...
                +0.5*r*(feval(f,h*i)+feval(f,h*(i-2))); %第二层
end
%显式迭代计算其他时间层
for j=3:m
    for i=2:n-1
        U(i,j)=2*(1-r)*U(i,j-1)+r*(U(i-1,j-1)+U(i+1,j-1))-U(i,j-2);
    end
end
U=U';
