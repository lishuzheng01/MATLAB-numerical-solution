%�Ľ��ǵ�˹Ԥ��-У������ʽMATLAB ����
function [t, y]=m4adams(df,tspan,y0,h)
%��;: �Ľ��ǵ�˹Ԥ��-У����ʽ�ⳣ΢�ַ���y'=f(t, y), y(t0)=y0
%��ʽ: [t, y]=m4adams(df,tspan,y0,h), dfΪ����f(t,y)���ʽ, 
%tspanΪ�������[t0,tn],y0Ϊ��ֵ, hΪ����, tΪ�ڵ�, yΪ��ֵ��
t=tspan(1) : h : tspan(2);
[~,y]=m4rkutta(df,[t(1),t(4)],y0,h); %���Ľ�����-������ʽ����4����ֵ
for n=4:(length(t)-1)
  p=y(n)+h/24*(55*feval(df,t(n),y(n))-59*feval(df,t(n-1),y(n-1)) ...
      +37*feval(df,t(n-2),y(n-2))-9*feval(df,t(n-3),y(n-3)));
  y(n+1)=y(n)+h/24*(feval(df,t(n-2),y(n-2))-5*feval(df,t(n-1),y(n-1))...
      +19*feval(df,t(n),y(n))+9*feval(df,t(n+1),p));
end