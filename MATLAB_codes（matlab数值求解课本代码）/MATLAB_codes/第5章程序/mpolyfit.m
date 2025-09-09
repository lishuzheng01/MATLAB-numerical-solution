%多项式拟合MATLAB程序
function p=mpolyfit(x,y,m)
%用途: 多项式拟合
%格式: p=mpolyfit(x,y,m), x,y为数据向量,m为拟合
%多项式的次数,p返回多项式系数降幂排列
A=zeros(m+1,m+1); b=zeros(m+1,1);
for i=0:m
    for j=0:m
        A(i+1,j+1)=sum(x.^(i+j));
    end
    b(i+1)=sum(x.^i.*y);
end
a=A\b; %解法方程组
p=fliplr(a');  %按降幂排列