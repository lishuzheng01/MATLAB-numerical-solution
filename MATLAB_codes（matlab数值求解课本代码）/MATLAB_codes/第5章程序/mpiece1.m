%�ֶ����Բ�ֵMATLAB����
function yy=mpiece1(x,y,xx)
%��;���ֶ����Բ�ֵ
%��ʽ��yy=mpiece1(x,y,xx), x�ǽڵ�����, y�ǽڵ��Ӧ�ĺ�
%��ֵ����, xx�ǲ�ֵ��(�����Ƕ��), yy���ز�ֵ���
n=length(x);
for j=1:length(xx)
    for i=2:n
        if xx(j)<x(i)
            yy(j)=y(i-1)*(xx(j)-x(i))/(x(i-1)-x(i)) ...
                     +y(i)*(xx(j)-x(i-1))/(x(i)-x(i-1));
            break;
        end
    end
end