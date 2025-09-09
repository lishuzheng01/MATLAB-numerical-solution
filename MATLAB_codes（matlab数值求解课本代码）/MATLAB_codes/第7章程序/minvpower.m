%ԭ��λ�Ʒ��ݷ�MATLAB����
function [lam,v,k]=minvpower(A,x,alpha,epsi,N)
%��;: ��ԭ��λ�Ʒ��ݷ�������ģ�������ֵ�Ͷ�Ӧ����������
%��ʽ: [lam,v,k]=minvpower(A,x,alpha,eps,N)
%����: AΪn�׷���,xΪ��ʼ����,epsiΪ���ƾ���(Ĭ��1.e-6),
%        NΪ����������(Ĭ��500��),alphaΪģ���Ľ�������ֵ.
%���: lam--����alpha���������ֵ,v--��Ӧ����������,k--��������
if nargin<5, N=500; end
if nargin<4,epsi=1e-6;end
m=0.5; k=0; err=1;
A=A-alpha*eye(length(x));
[L,U,P]=lu(A);
while (k<N) & (err>epsi)
   [~,t]=max(abs(x));
   m1=x(t); v=x/m1;
   z=L\(P*v);  x=U\z;
   err=abs(1/m1-1/m);
   k=k+1;  m=m1;
end
lam=alpha+1/m;