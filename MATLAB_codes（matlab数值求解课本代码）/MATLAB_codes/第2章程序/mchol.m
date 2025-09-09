%Cholesky�ֽⷨMATLAB����
function [x,l,d]=mchol(A,b)
%��;: ��Cholesky�ֽⷨ�ⷽ����Ax=b
%��ʽ: [x,l,d]=mchol(A,b), AΪϵ������, bΪ�Ҷ�����
%����: x--������, l--��λ��������, d--����D�ĶԽ��� 
format short
%LDL'�ֽ�
n=length(b); d=zeros(1,n); 
l=eye(n,n);  t=zeros(n,n);
d(1)=A(1,1); t(1,1)=A(1,1);
for i=2 : n
    t(i,1)=A(i,1); l(i,1)=t(i,1)/d(1);
end
for k=2 : n
    s=0;
    for r=1:k-1
        s=s+t(k,r)*l(k,r);
    end
    d(k)=A(k,k)-s;
    for i=k+1 : n
        s=0;
        for r=1:k-1
            s=s+t(i,r)*l(k,r);
        end
        t(i,k)=A(i,k)-s;
        l(i,k)=t(i,k)/d(k);
    end
end
%��������Ƿ����� Ly=b 
y=zeros(n,1);
y(1)=b(1); 
for i=2:n
    y(i)=b(i)-l(i,1:i-1)*y([1:i-1]);
end
%���ԽǷ����� Dz=y 
for i=1:n
    z(i)=y(i)/d(i);
end
%��������Ƿ����� L'x=z
u=l'; x=zeros(n,1);
x(n)=z(n);
for i=(n-1):-1:1
    x(i)=z(i)-u(i,i+1:n)*x(i+1:n);
end
