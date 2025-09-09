%�����������MATLAB����
function [x,iter]=mfpm(phi,x0,ep,N)
%��;: �ò����������������Է���f(x)=0�и�����[a,b]�е�һ����
%��ʽ: [x,iter]=mfpm(phi,x0,ep,N), phiΪ��������, x0Ϊ��ֵ, epΪ����
%(Ĭ��1e-5), NΪ����������(Ĭ��500), x,iter�ֱ�Ϊ���Ƹ��͵�������
if nargin<4 N=500;end
if nargin<3 ep=1e-6;end
iter=0;
while iter<N
   x=feval(phi,x0);
   if abs(x-x0)<ep
      break;
   end
   x0=x; iter=iter+1;
end