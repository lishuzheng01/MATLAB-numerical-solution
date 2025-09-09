%拉格朗日插值法MATLAB程序
function yp=mlagr(x,y,xp)
%用途：拉格朗日插值法
%格式：yp=mlagr(x,y,xp), x是节点向量, y是节点对应的函
%数值向量, xp是插值点(可以是多个), yp返回插值结果
n=length(x); m=length(xp);
yp=zeros(1,m);  c1=ones(n-1,1); c2=ones(1,m);
for i=1:n
   xb=x([1:i-1,i+1:n]);
   yp=yp+y(i)*prod((c1*xp-xb'*c2)./(x(i)-xb'*c2));
end