%Householder�任MATLAB����
function H=mhouseh(x)
%��;: ��������x,����Householder�任H,ʹ��Hx=(*,0,...,0)'
%��ʽ: H=mhouseh(x)
%xΪ����������, H����Householder�任����
n=length(x); I=eye(n); sn=sign(x(1));
if sn==0, sn=1; end
z=x(2:n);
if(norm(z,inf)==0), H=I; return; end
a=sn*norm(x,2);
u=x; u(1)=u(1)+a;
gamma=1.0/(a*(a+x(1)));
H=I-gamma*u*u';