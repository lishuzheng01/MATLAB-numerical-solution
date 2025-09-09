%����ʽ���MATLAB����
function p=mpolyfit(x,y,m)
%��;: ����ʽ���
%��ʽ: p=mpolyfit(x,y,m), x,yΪ��������,mΪ���
%����ʽ�Ĵ���,p���ض���ʽϵ����������
A=zeros(m+1,m+1); b=zeros(m+1,1);
for i=0:m
    for j=0:m
        A(i+1,j+1)=sum(x.^(i+j));
    end
    b(i+1)=sum(x.^i.*y);
end
a=A\b; %�ⷨ������
p=fliplr(a');  %����������