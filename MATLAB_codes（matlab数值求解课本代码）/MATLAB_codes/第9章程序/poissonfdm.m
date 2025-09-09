%程序: 泊松方程有限差分法程序 
function pois_error=poissonfdm(n)
%五点格式计算泊松方程的数值解
%-(u_{xx}+u_{yy}) = f in [0,1]x[0,1]
%f =2*sin(pi*x)*sin(pi*y); Dirichlet边界条件
%解析解: u=(1/pi^2)*sin(pi*x)*sin(pi*y)
clc
x=1/n*(0:n);  %x方向坐标
y=1/n*(0:n);  %y方向坐标
%下边界go: y=0
go=zeros(1,n+1);  %u0=0
%左边界gl: x=0
gl=zeros(1,n+1);  %u0=(1/pi^2)*sin(pi*x)*sin(pi*y)
%右边界gr: x=1
gr=zeros(1,n+1);
%上边界gu: y=1
gu=zeros(1,n+1);
%**************计算向量f和p ****************
f=[ ];
x=1/n*(1:n-1); %内节点横坐标
y=1/n*(1:n-1); %内节点纵坐标
a=2*sin(pi*x);
for j=1:n-1   
    f =[f,a*sin(pi*y(j))];
end
%p是边界处理后的线性方程组A*u=p的右端向量
p=f;
%next to the lower border
p(1:n-1) = p(1:n-1)  + n^2*go(2:n); 
%next to the left border
p(1:n-1:(n-1)*(n-2)+1)  = p(1:n-1:(n-1)*(n-2)+1)  + n^2*gl(2:n);
%next to the right border
p(n-1:n-1:(n-1)^2) = p(n-1:n-1:(n-1)^2)   + n^2*gr(2:n);
%next to the upper border
p((n-1)*(n-2)+1:(n-1)^2)=p((n-1)*(n-2)+1:(n-1)^2)+n^2*gu(2:n);
%
p=p';
%***************计算矩阵A **********************
e=ones((n-1)^2,1);
A=spdiags([-e -e 4*e -e -e],[-(n-1) -1 0 1 n-1],(n-1)^2,(n-1)^2);
for j=1:n-2   
   A(j*(n-1),j*(n-1)+1)=0;   A(j*(n-1)+1,j*(n-1))=0;
end
A=A*n^2;
%***************计算解向量u **************************
u=A\p;  u=u';
%***************比较数值解与解析解****************
anal=[]; %计算解析解
a=(1/pi^2)*sin(pi*x);
for j=1:n-1   
    anal=[anal,a*sin(pi*y(j))];
end
pois_error=max(abs(anal-u)); %计算最大偏差
%***************可视化数值解和解析解 **********************
% order the solution vectors as matrices
sol=reshape(anal,n-1,n-1); %解析解向量重排成矩阵
sol=sol';
mat=reshape(u,n-1,n-1); %数值解向量u重排为矩阵
mat=mat';
figure
subplot(2,1,1)
surf([0:1/n:1],[0:1/n:1],[go;gl(2:n)',mat,gr(2:n)';gu])
axis([0 1 0 1 0 0.1]);%view(3);      
%shading interp
xlabel('x-direction'), ylabel('y-directon')
title('数值解'), colorbar
subplot(2,1,2)
surf([0:1/n:1],[0:1/n:1],[go;gl(2:n)',sol,gr(2:n)';gu])
axis([0 1 0 1 0 0.1]);%view(3);            
%shading interp
xlabel('x-direction'), ylabel('y-directon')
title('解析解'),colorbar
