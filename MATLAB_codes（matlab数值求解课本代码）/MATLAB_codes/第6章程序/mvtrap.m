%�䲽�����ι�ʽ~MATLAB ����
function [T,k,n]=mvtrap(f,a,b,epsi)
%��;:�ñ䲽�����ι�ʽ�����.
%��ʽ:[T,k,n]=mvtrap(f,a,b,epsi),fΪ��������,a,bΪ���������
%���Ҷ˵�,epsiΪ���ƾ���,T���ػ��ֽ���ֵ,kΪ�Է�����Ĵ���,
%nΪ����ĵȷ���
h=b-a;
T1=h*(feval(f,a)+feval(f,b))/2;
T=T1/2+(h/2)*feval(f,a+h/2);
k=1;n=2;
while(abs(T-T1)>epsi)
    h=h/2; T1=T;
    T=T1/2+(h/2)*sum(feval(f,a+h/2:h:b-h/2));
    k=k+1; n=2*n;
end