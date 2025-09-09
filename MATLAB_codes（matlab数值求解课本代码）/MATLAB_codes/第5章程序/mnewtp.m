%牛顿插值法MATLAB程序
function yy=mnewtp(x,y,xx)
%用途：牛顿插值法
%格式：yy=mnewtp(x,y,xx), x是节点向量, 
%y是节点对应的函数值向量, xx是插值点, yy返回插值结果
n=length(x);  yy=y(1);
y1=0; lx=1;
for i=1:n-1
   for j=i+1:n
      y1(j)=(y(j-1)-y(j))/(x(j-i)-x(j));  %计算差商
   end
   c(i)=y1(i+1); lx=lx*(xx-x(i));
   yy=yy+c(i)*lx;  %计算牛顿插值多项式的值
   y=y1;
end
