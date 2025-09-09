%Jacobi����MATLAB����
function [D,V]=mejacobi(A,epsi)
%��;: ��Jacobi������ʵ�Գƾ���A������ֵ����������
%��ʽ: [D,V]=mejacobi(A,epsi),
%����: AΪn�׶ԳƷ���, epsiΪ�������(Ĭ��1.e-6) 
%���: D�ǶԽ���, ��Խ�ԪΪA��n������ֵ, V����������������
if nargin<2,epsi=1e-6;end
n=size(A,1); D=zeros(n); V=eye(n);
%����A�ķǶԽ�Ԫ����ֵ���Ԫ�����ڵ���p����q
[w1,p]=max(abs(A-diag(diag(A))));
[~,q]=max(w1); p=p(q);
while (1)
    d=(A(q,q)-A(p,p))/(2*A(p,q));
    if(d>0)
        t=-d+sqrt(d^2+1);
    else if(d<0)
            t=-d-sqrt(d^2+1);
        else
            t=0;
        end
    end
    c=1/sqrt(t^2+1);  s=c*t;
    G=[c s; -s c];
    A([p q],:)=G'*A([p q],:);
    A(:,[p q])=A(:,[p q])*G;
    V(:,[p q])=V(:,[p q])*G;
    [w1,p]=max(abs(A-diag(diag(A))));
    [~,q]=max(w1);  p=p(q);
    if (abs(A(p,q))<epsi*sqrt(sum(diag(A).^2)/n))
        break;
    end
end
D=diag(diag(A));