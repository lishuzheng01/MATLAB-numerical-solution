%高阶插值的Runge现象
xx=-5:0.05:5;y=1./(1+xx.^2);
x1=-5:2:5;y1=1./(1+x1.^2);
x2=-5:1:5; y2=1./(1+x2.^2);
yy1=mlagr(x1,y1,xx);
yy2=mlagr(x2,y2,xx);
plot(xx,yy1,'k--'); hold on
plot(xx,yy2,'kx'); plot(xx,y,'k');
legend('5次多项式插值','10次多项式插值','被插函数的图形');
axis([-5 5 -0.5 2])