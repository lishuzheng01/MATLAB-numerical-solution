%�ݶ��½���MATLAB����
function [x,iter]=mgrad(A,b,x,ep,N)
%��;�����ݶ��½��������Է�����Ax=b
%��ʽ��[x,iter]=mgrad(A,b,x,ep,N)  ����AΪϵ������, bΪ�Ҷ���
%��, xΪ��ʼ����(Ĭ��������), epΪ����(Ĭ��1e-6),NΪ��������
%��(Ĭ��1000��), ���ز���x,iter�ֱ�Ϊ���ƽ������͵�������
if nargin<5, N=1000; end
if nargin<4, ep=1e-6; end
if nargin<3, x=zeros(size(b)); end
for iter=1:N
   r=b-A*x;
   if norm(r)<ep, break; end
   a=r'*r/(r'*A*r);
   x=x+a*r;
end