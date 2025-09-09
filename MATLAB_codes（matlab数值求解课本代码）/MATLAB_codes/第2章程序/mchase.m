%追赶法MATLAB程序
function x=mchase(a,b,c,d)
%用途：追赶法解三对角方程组
%格式：x= mchase(a,b,c,d),  a为次下对角线元素向量，b主对角元素
%向量，c为次上对角线元素向量，d为右端向量，x返回解向量
n=length(b);
for k=2:n
    b(k)=b(k)-a(k)/b(k-1)*c(k-1);
    d(k)=d(k)-a(k)/b(k-1)*d(k-1);
end
x(n)=d(n)/b(n);
for k=n-1:-1:1
    x(k)=(d(k)-c(k)*x(k+1))/b(k);
end
