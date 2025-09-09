%SOR������MATLAB����
function [x,iter]=msor(A,b,omega,x,ep,N)
%��;:��SOR�����������Է�����Ax=b
%��ʽ:[x,iter]=msor(A,b,omega,x,ep,N)  AΪϵ������, bΪ��
%������, omegaΪ�ɳ�����(Ĭ��1.05), xΪ��ʼ����(Ĭ��������),
%epΪ����(Ĭ��e-6), NΪ����������(Ĭ��500��), ���ز���x,
%iter�ֱ�Ϊ���ƽ�͵�������
if nargin<6, N=500; end
if nargin<5, ep=1e-6; end
if nargin<4, x=zeros(size(b)); end
if nargin<3, omega=1.05; end
D=diag(diag(A)); L=D-tril(A); U=D-triu(A);
for iter=1:N
   x=(D-omega*L)\(((1-omega)*D+omega*U)*x+omega*b);
   err=norm(b-A*x)/norm(b);
   if err<ep, break; end
end