%分段线性插值MATLAB程序
function yy=mpiece1(x,y,xx)
%用途：分段线性插值
%格式：yy=mpiece1(x,y,xx), x是节点向量, y是节点对应的函
%数值向量, xx是插值点(可以是多个), yy返回插值结果
n=length(x);
for j=1:length(xx)
    for i=2:n
        if xx(j)<x(i)
            yy(j)=y(i-1)*(xx(j)-x(i))/(x(i-1)-x(i)) ...
                     +y(i)*(xx(j)-x(i-1))/(x(i)-x(i-1));
            break;
        end
    end
end