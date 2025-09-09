%main_heat program
f=@(x)4*x-4*x.^2;
c1=0; c2=0; c=1; a=1; b=0.2; 
n=6;m=11;
U=forwarddif(f,c1,c2,c,a,b,n,m)
surf(U);
xlabel('x-Öá'), ylabel('t-Öá')
