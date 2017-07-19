syms t
t=-10:0.1:10; 
b=(1/sqrt(2*pi)).*exp(-t.*t/2);
plot(t,b)
title('高斯函数');
         
Db1 = diff(b)     %一阶导数
Db2 = diff(b,2)   %二阶求导
 
figure
plot(Db1)
title('一阶导数');
 
figure
plot(Db2)
title('二阶导数');
