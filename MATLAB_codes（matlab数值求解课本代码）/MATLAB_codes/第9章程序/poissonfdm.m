%����: ���ɷ������޲�ַ����� 
function pois_error=poissonfdm(n)
%����ʽ���㲴�ɷ��̵���ֵ��
%-(u_{xx}+u_{yy}) = f in [0,1]x[0,1]
%f =2*sin(pi*x)*sin(pi*y); Dirichlet�߽�����
%������: u=(1/pi^2)*sin(pi*x)*sin(pi*y)
clc
x=1/n*(0:n);  %x��������
y=1/n*(0:n);  %y��������
%�±߽�go: y=0
go=zeros(1,n+1);  %u0=0
%��߽�gl: x=0
gl=zeros(1,n+1);  %u0=(1/pi^2)*sin(pi*x)*sin(pi*y)
%�ұ߽�gr: x=1
gr=zeros(1,n+1);
%�ϱ߽�gu: y=1
gu=zeros(1,n+1);
%**************��������f��p ****************
f=[ ];
x=1/n*(1:n-1); %�ڽڵ������
y=1/n*(1:n-1); %�ڽڵ�������
a=2*sin(pi*x);
for j=1:n-1   
    f =[f,a*sin(pi*y(j))];
end
%p�Ǳ߽紦�������Է�����A*u=p���Ҷ�����
p=f;
%next to the lower border
p(1:n-1) = p(1:n-1)  + n^2*go(2:n); 
%next to the left border
p(1:n-1:(n-1)*(n-2)+1)  = p(1:n-1:(n-1)*(n-2)+1)  + n^2*gl(2:n);
%next to the right border
p(n-1:n-1:(n-1)^2) = p(n-1:n-1:(n-1)^2)   + n^2*gr(2:n);
%next to the upper border
p((n-1)*(n-2)+1:(n-1)^2)=p((n-1)*(n-2)+1:(n-1)^2)+n^2*gu(2:n);
%
p=p';
%***************�������A **********************
e=ones((n-1)^2,1);
A=spdiags([-e -e 4*e -e -e],[-(n-1) -1 0 1 n-1],(n-1)^2,(n-1)^2);
for j=1:n-2   
   A(j*(n-1),j*(n-1)+1)=0;   A(j*(n-1)+1,j*(n-1))=0;
end
A=A*n^2;
%***************���������u **************************
u=A\p;  u=u';
%***************�Ƚ���ֵ���������****************
anal=[]; %���������
a=(1/pi^2)*sin(pi*x);
for j=1:n-1   
    anal=[anal,a*sin(pi*y(j))];
end
pois_error=max(abs(anal-u)); %�������ƫ��
%***************���ӻ���ֵ��ͽ����� **********************
% order the solution vectors as matrices
sol=reshape(anal,n-1,n-1); %�������������ųɾ���
sol=sol';
mat=reshape(u,n-1,n-1); %��ֵ������u����Ϊ����
mat=mat';
figure
subplot(2,1,1)
surf([0:1/n:1],[0:1/n:1],[go;gl(2:n)',mat,gr(2:n)';gu])
axis([0 1 0 1 0 0.1]);%view(3);      
%shading interp
xlabel('x-direction'), ylabel('y-directon')
title('��ֵ��'), colorbar
subplot(2,1,2)
surf([0:1/n:1],[0:1/n:1],[go;gl(2:n)',sol,gr(2:n)';gu])
axis([0 1 0 1 0 0.1]);%view(3);            
%shading interp
xlabel('x-direction'), ylabel('y-directon')
title('������'),colorbar
