%��������ɭ��ʽMATLAB����
function s=msimp(f,a,b,n)
%��;: �ø�������ɭ��ʽ�����.
%��ʽ: s=msimp(f,a,b,n)  fΪ��������;a,bΪ������
%������Ҷ˵�,nΪ����ĵȷ���,s���ػ��ֽ���ֵ
h=(b-a)/n;
x=linspace(a,b,2*n+1);
y=feval(f,x);
s1=sum(y(3:2:2*n-1));
s2=sum(y(2:2:2*n));
s=(h/6)*(y(1)+2*s1+4*s2+y(n+1));