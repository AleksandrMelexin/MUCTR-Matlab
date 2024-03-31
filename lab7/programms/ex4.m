x=[0,1.13,1.5,2.25,3];
y=[4.57,0.68,0.39,-1.9,-4.4];
% Вычислим приближения с различной степенью
p0=polyfit(x,y,0);
p1=polyfit(x,y,1);
p2=polyfit(x,y,2);
p3=polyfit(x,y,3);
% Вычислим ошибки (СКО)
y0=polyval(p0,x);
y1=polyval(p1,x);
y2=polyval(p2,x);
y3=polyval(p3,x);
err0=(1/5*sum((y-y0).^2))^0.5
err1=(1/5*sum((y - y1).^2))^0.5
err2=(1/5*sum((y-y2).^2))^0.5
err3=(1/5*sum((y-y3).^2))^0.5
