%�����е㹫ʽMATLAB����
function s=midpt(f,a,b,n)
%��;���ø����е㹫ʽ�����.
%��ʽ��s=mintm(f,a,b,n),fΪ��������,a,bΪ������
%������Ҷ˵�,nΪ����ĵȷ���,s���ػ��ֽ���ֵ
h=(b-a)/n;
x=linspace(a+h/2,b-h/2,n);
y=feval(f,x);
s=h*sum(y);