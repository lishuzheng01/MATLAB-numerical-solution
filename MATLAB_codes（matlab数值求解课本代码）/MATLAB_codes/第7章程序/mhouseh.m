%Householder变换MATLAB程序
function H=mhouseh(x)
%用途: 对于向量x,构造Householder变换H,使得Hx=(*,0,...,0)'
%格式: H=mhouseh(x)
%x为输入列向量, H返回Householder变换矩阵
n=length(x); I=eye(n); sn=sign(x(1));
if sn==0, sn=1; end
z=x(2:n);
if(norm(z,inf)==0), H=I; return; end
a=sn*norm(x,2);
u=x; u(1)=u(1)+a;
gamma=1.0/(a*(a+x(1)));
H=I-gamma*u*u';