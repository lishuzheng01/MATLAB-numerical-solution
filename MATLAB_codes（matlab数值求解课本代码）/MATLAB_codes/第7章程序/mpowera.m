%ԭ��λ���ݷ�����MATLAB����
function [lam,v,k]=mpowera(A,x,alpha,epsi,N)
%��;: ��ԭ��λ���ݷ�������ģ�������ֵ�Ͷ�Ӧ����������
%��ʽ: [lam,v,k]=mpowera(A,x,alpha,epsi,N)
%����: AΪn�׷���, xΪ��ʼ����, epsiΪ���ƾ���(Ĭ��10e-6), 
%N����������(Ĭ��500��), alphaΪԭ��λ�Ʋ���. 
%���: lam--��ģ��������ֵ, v--��Ӧ����������, k--��������.
if nargin<5,N=500;end
if nargin<4,epsi=1e-6;end
m=0; k=0; err=1;
A=A-alpha*eye(length(x));
while (k<N) & (err>epsi)
    v=A*x; [~,t]=max(abs(v));
    m1=v(t); x=v/m1;
    err=abs(m1-m);
    m=m1; k=k+1;
end
lam=m1+alpha; v=x;