%����QR����MATLAB���� 
function [Iter,D]=mqrdec(A,epsi)
%��;: �û���QR�㷨��ʵ�����ȫ������ֵ.
%����: n��ʵ����A, ���ƾ���epsi(Ĭ����1.e-5)
%���: ��������Iter, A��ȫ������ֵD
%���ú���: mhessenb.m, mhessenqr.m, eig--������1,2����
if nargin<2,epsi=1e-5;end
n=size(A,1); D=zeros(n,1); 
i=n; Iter=0;  %��ʼ��
A=mhessenb(A);  %������AΪHessenberg��
%�û���QR�㷨���е���
while (1)
    if n<=2
        la=eig(A(1:n,1:n)); D(1:n)=la';  break;
    end
    Iter=Iter+1;
    [Q,R]=mhessenqr(A); %����Hessenberg ������QR�ֽ�
    A=R*Q;  %���������Ʊ任
    %����ĳ�������ж��Ƿ���ֹ
    for k=n-1:-1:1
        if abs(A(k+1,k))<epsi
            if n-k<=2
                la=eig(A(k+1:n,k+1:n));
                j=i-n+k+1; D(j:i)=la';
                i=j-1; n=k; break;
            end
        end
    end
end