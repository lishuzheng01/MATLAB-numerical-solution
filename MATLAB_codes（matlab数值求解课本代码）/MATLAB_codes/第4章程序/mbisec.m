%���ַ�MATLAB����
function [x,iter]=mbisec(f,a,b,ep)
%��;: �ö��ַ�������Է���f(x)=0�и�����[a,b]�е�һ����
%��ʽ: [x,iter]=mbisec(f,a,b,ep),  fΪ�������ʽ��a,bΪ����
%���Ҷ˵�, epΪ����, x, iter�ֱ𷵻ؽ��Ƹ��Ͷ��ִ���
x=(a+b)/2.0;  iter=0;
while abs(feval(f,x))>ep|(b-a>ep)
    if feval(f,x)*feval(f,a)<0
        b=x;
    else
        a=x;
    end
    x=(a+b)/2.0; 
    iter=iter+1;
end