%���ؽ�-ʷ����ɭ���ٷ�MATLAB����
function [x,iter]=maitken(phi,x0,ep,N)
%��;: �ð��ؽ�-ʷ����ɭ���ٷ���f(x)=0�Ľ�
%��ʽ: [x,iter]=maitken(phi,x0,ep,N), phiΪ��������, x0Ϊ������ֵ,epΪ��
%��(Ĭ��1e-5), NΪ����������(Ĭ��500), x,iter�ֱ𷵻ؽ��Ƹ��͵�������
if nargin<4,N=500;end
if nargin<3,ep=1e-5;end
iter=0;
while iter<N
   y=feval(phi,x0);
   z=feval(phi,y);
   x=x0-(y-x0)^2/(z-2*y+x0);
   if abs(x-x0)<ep, break;  end
   x0=x; iter=iter+1;
end