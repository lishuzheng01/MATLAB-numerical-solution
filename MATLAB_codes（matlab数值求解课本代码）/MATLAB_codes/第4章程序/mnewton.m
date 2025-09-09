%ţ�ٷ�MATLAB����
function [x,iter]=mnewton(f,df,x0,ep,N)
%��;: ��ţ�ٷ��������Է���f(x)=0
%��ʽ: [x,iter]=mnewton(f,df,x0,ep,N)
%f��df�ֱ�Ϊ��ʾf(x)���䵼��, x0Ϊ������ֵ, epΪ����(Ĭ��1e-5), 
%NΪ����������(Ĭ��Ϊ500), x,iter�ֱ𷵻ؽ��Ƹ��͵�������
if nargin<5,N=500;end
if nargin<4,ep=1e-5;end
iter=0;
while iter<N
    x=x0-feval(f,x0)/feval(df,x0);
    if abs(x-x0)<ep
        break;
    end
    x0=x; iter=iter+1;
end