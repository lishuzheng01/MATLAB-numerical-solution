%���߷�MATLAB����
function [x,iter]=msecant(f,x0,x1,ep,N)
%��;:�ø��߷��������Է���f(x)=0
%��ʽ:[x,iter]=msecant(f,x0,x1,ep,N)  fΪf(x)�ı��ʽ, x0, x1Ϊ
%������ֵ, epΪ����(Ĭ��1e-5), NΪ����������(Ĭ��Ϊ500), 
%x,iter�ֱ𷵻ؽ��Ƹ��͵�������
if nargin<5,N=500;end
if nargin<4,ep=1e-5;end
iter=0;
while iter<N
    x=x1-(x1-x0)*feval(f,x1)/(feval(f,x1)-feval(f,x0));
    if abs(x-x1)<ep,  break;  end
    x0=x1; x1=x; iter=iter+1;
end