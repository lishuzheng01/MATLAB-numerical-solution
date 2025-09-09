%Jacobi方法MATLAB程序
function [D,V]=mejacobi(A,epsi)
%用途: 用Jacobi方法求实对称矩阵A的特征值和特征向量
%格式: [D,V]=mejacobi(A,epsi),
%输入: A为n阶对称方阵, epsi为容许误差(默认1.e-6) 
%输出: D是对角阵, 其对角元为A的n个特征值, V返回特征向量矩阵
if nargin<2,epsi=1e-6;end
n=size(A,1); D=zeros(n); V=eye(n);
%计算A的非对角元绝对值最大元素所在的行p和列q
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