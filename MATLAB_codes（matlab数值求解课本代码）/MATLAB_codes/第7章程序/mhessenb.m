%������Hessenberg����MATLAB����
function A=mhessenb(A)
%��;: ��Householder�任������Ϊ��Hessenberg����.
%����: n��ʵ����A; ���: A����Hessenberg��
%���ú���: mhouseh.m
n = size(A,1);
for k=1:(n-2)
    x=A(k+1:n,k); H=mhouseh(x);
    A(k+1:n,1:n)=H*A(k+1:n,1:n);
    A(1:n,k+1:n)=A(1:n,k+1:n)*H;
end