%ţ�ٲ�ֵ��MATLAB����
function yy=mnewtp(x,y,xx)
%��;��ţ�ٲ�ֵ��
%��ʽ��yy=mnewtp(x,y,xx), x�ǽڵ�����, 
%y�ǽڵ��Ӧ�ĺ���ֵ����, xx�ǲ�ֵ��, yy���ز�ֵ���
n=length(x);  yy=y(1);
y1=0; lx=1;
for i=1:n-1
   for j=i+1:n
      y1(j)=(y(j-1)-y(j))/(x(j-i)-x(j));  %�������
   end
   c(i)=y1(i+1); lx=lx*(xx-x(i));
   yy=yy+c(i)*lx;  %����ţ�ٲ�ֵ����ʽ��ֵ
   y=y1;
end
