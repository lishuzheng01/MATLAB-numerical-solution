%�������ι�ʽMATLAB����
function s=mtrap(f,a,b,n)
%��;: �ø������ι�ʽ�����.
%��ʽ: s=mtrap(f,a,b,n),  fΪ��������,a,bΪ������
%������Ҷ˵�,nΪ����ĵȷ���,s���ػ��ֽ���ֵ
h=(b-a)/n;
x=linspace(a,b,n+1);
y=feval(f,x);
s=0.5*h*(y(1)+2*sum(y(2:n))+y(n+1));