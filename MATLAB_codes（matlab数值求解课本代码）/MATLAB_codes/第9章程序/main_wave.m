%main_wave program
clear all
f=@(x)sin(pi*x)+sin(2*pi*x); g=@(x)0;
c=2; a=1; b=0.5; 
n=21;m=21;
U=finitedif(f,g,c,a,b,n,m); %������ֵ��
figure (1)
surf(U);  title('��ֵ��'), 
xlabel('x-��'), ylabel('t-��')
h=a/(n-1); k=b/(m-1);
for i=1:n,  x(i)=i*h; end
for j=1:m, t(j)=j*k; end
for j=1:m  %���������
    for i=1:n
        U1(i,j)=sin(pi*x(i))*cos(2*pi*t(j))+sin(2*pi*x(i))*cos(4*pi*t(j));
    end
end
U1=U1'; 
figure (2)
surf(U1); title('������')
xlabel('x-��'), ylabel('t-��')
err=norm(U-U1,'fro')
