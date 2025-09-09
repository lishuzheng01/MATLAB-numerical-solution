%�ݷ�MATLAB����
function [lam,v,k]=mepower(A,x,epsi,N)
%��;: ���ݷ�������ģ�������ֵ�Ͷ�Ӧ����������
%��ʽ: [lam,v,k]=mepower(A,x,epsi,N)
%����: AΪn�׷���, xΪ��ʼ����, epsi���ƾ���, NΪ����������.
%���: lam--��ģ��������ֵ, v--��Ӧ����������, k--��������.
if nargin<4, N=500; end
if nargin<3, epsi=1e-6; end
m=0; k=0; err=1;
while (k<N) & (err>epsi)
   v=A*x;
   [~,t]=max(abs(v));
   m1=v(t);  x=v/m1;
   err=abs(m1-m);
   m=m1;   k=k+1;
end
lam=m1; v=x;