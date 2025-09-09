function df=mdf(t,y)
alpha=30; beta=2.8; sigma=12;
df=zeros(3,1);
df(1)=-sigma*y(1)+sigma*y(2);
df(2)=alpha*y(1)-y(2)-y(1)*y(3);
df(3)=y(1)*y(2)-beta*y(3);
end