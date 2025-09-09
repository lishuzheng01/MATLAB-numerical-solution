%�������ղ�ֵ��MATLAB����
function yp=mlagr(x,y,xp)
%��;���������ղ�ֵ��
%��ʽ��yp=mlagr(x,y,xp), x�ǽڵ�����, y�ǽڵ��Ӧ�ĺ�
%��ֵ����, xp�ǲ�ֵ��(�����Ƕ��), yp���ز�ֵ���
n=length(x); m=length(xp);
yp=zeros(1,m);  c1=ones(n-1,1); c2=ones(1,m);
for i=1:n
   xb=x([1:i-1,i+1:n]);
   yp=yp+y(i)*prod((c1*xp-xb'*c2)./(x(i)-xb'*c2));
end