%��˹-���¶�������MATLAB����
function [x,iter]=mseidel(A,b,x,ep,N)
%��;���ø�˹-���¶������������Է�����Ax=b
%��ʽ��[x,iter]=mseidel(A,b,x,ep,N)  AΪϵ������, bΪ�Ҷ���
%��, xΪ��ʼ����(Ĭ��������), epΪ����(Ĭ��1e-6), NΪ����
%������(Ĭ��500��), ���ز���x,iter�ֱ�Ϊ���ƽ������͵�������
if nargin<5, N=500; end
if nargin<4, ep=1e-6; end
if nargin<3, x=zeros(size(b)); end
D=diag(diag(A)); L=D-tril(A); U=D-triu(A);
for iter=1:N
   x=(D-L)\(U*x+b);
   err=norm(b-A*x)/norm(b);
   if err<ep, break; end
end