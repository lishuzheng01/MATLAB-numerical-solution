%�����ݶȷ�MATLAB����
function [x,iter]=mcg(A,b,x,ep,N)
%��;���ù����ݶȷ������Է�����Ax=b
%��ʽ��[x,iter]=mcg(A,b,x,ep,N)  ���� AΪϵ������, bΪ�Ҷ�
%����, xΪ��ʼ����(Ĭ��������), epΪ����(Ĭ��1e-6),NΪ����
%������(Ĭ��500��), ���ز���x,iter�ֱ�Ϊ���ƽ������͵�������
if nargin<5, N=500; end
if nargin<4, ep=1e-6; end
if nargin<3, x=zeros(size(b)); end
r=b-A*x;
for iter=1:N
   rr=r'*r;
   if iter==1
       p=r;
   else
       beta=rr/rr1;
       p=r+beta*p;
   end
   q=A*p;
   alpha=rr/(p'*q);
   x=x+alpha*p;
   r=r-alpha*q;
   rr1=rr;
   if (norm(r)<ep), break; end
end