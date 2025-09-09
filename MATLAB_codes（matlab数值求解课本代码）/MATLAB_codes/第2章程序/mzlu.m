%����ԪLU�ֽⷨMATLAB����
function [x,L,U,P]=mzlu(A,b)
%��;: ������ԪLU�ֽⷨ�ⷽ����Ax=b
%��ʽ: [x,L,U,P]=mzlu(A,b),  AΪϵ������, bΪ�Ҷ�����,
%����: x-������, L-��λ�����Ǿ���, U-�����Ǿ���,
%P-ѡ��Ԫʱ��¼�н������û���
n=length(b); 
P=eye(n); %P��¼ѡ����Ԫʱ�������е��б任
%����ԪLU�ֽ�
for k=1:n
    A(k:n,k)=A(k:n,k)-A(k:n,1:k-1)*A(1:k-1,k);
    [s,m]=max(abs(A(k:n,k)));   %ѡ����Ԫ
    m=m+k-1;
    if m~=k
        A([k m],:)=A([m k],:);
        P([k m],:)=P([m k],:);
    end
    A(k+1:n,k)=A(k+1:n,k)/A(k,k);
    A(k,k+1:n)=A(k,k+1:n)-A(k,1:k-1)*A(1:k-1,k+1:n);
end
L=tril(A,-1)+eye(n,n); U=triu(A);
%�ⵥλ�����Ǿ��� Ly=b
newb=P*b;  y=zeros(n,1);
for k=1: n
    j=1: k-1;
    y(k)=newb(k)-L(k,j)*y(j);
end 
%�������Ƿ�����Ux=y
x=zeros(n,1);
for k=n:-1:1
     j=k+1:n;
     x(k)=(y(k)-U(k,j)*x(j))/U(k,k);
end