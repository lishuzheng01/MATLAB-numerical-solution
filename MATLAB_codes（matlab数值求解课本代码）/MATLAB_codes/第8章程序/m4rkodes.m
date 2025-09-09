%�Ľ׾������������ʽ��һ�׳�΢�ַ�����MATLAB����
function [t,y]=m4rkodes(df,tspan,y0,h)
%��;: �Ľ׾�������-������ʽ�ⳣ΢�ַ�����y'=f(t,y), y(t0)=y0
%��ʽ: [t,y]=m4rkodes(df,tspan,y0,h),dfΪ��������f(t,y)���ʽ, 
%tspanΪ�������[t0,tn],y0Ϊ��ֵ����, hΪ����, tΪ�ڵ�, yΪ������
t=tspan(1):h:tspan(2);
y=zeros(length(y0),length(t));
y(:,1)=y0(:);
 for n=1:(length(t)-1)
     k1=feval(df, t(n), y(:,n));
     k2=feval(df, t(n)+h/2, y(:,n)+h/2*k1);
     k3=feval(df, t(n)+h/2, y(:,n)+h/2*k2);
     k4=feval(df, t(n+1), y(:,n)+h*k3);
     y(:,n+1)=y(:,n)+h*(k1+2*k2+2*k3+k4)/6;
 end