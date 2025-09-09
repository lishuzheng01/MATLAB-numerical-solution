%�Ľ���ŷ������MATLAB����
function [t,y]=meuler(df,tspan,y0,h)
%��;: �Ľ���ŷ�������ⳣ΢�ַ���y'=f(t,y), y(t0)=y0
%��ʽ: [t,y]=meuler(df,tspan,y0,h), dfΪ����f(t,y), tspanΪ���
%����[t0,tn], y0Ϊ��ֵy(t0), hΪ����, [t,y]���ؽڵ����ֵ�����
t=tspan(1):h:tspan(2);  y(1)=y0;
for n=1:(length(t)-1)
    k1=feval(df,t(n),y(n));
    y(n+1)=y(n)+h*k1;
    k2=feval(df,t(n+1),y(n+1));
    y(n+1)=y(n)+h*(k1+k2)/2;
end