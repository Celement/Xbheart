syms t
t=-10:0.1:10; 
b=(1/sqrt(2*pi)).*exp(-t.*t/2);
plot(t,b)
title('��˹����');
         
Db1 = diff(b)     %һ�׵���
Db2 = diff(b,2)   %������
 
figure
plot(Db1)
title('һ�׵���');
 
figure
plot(Db2)
title('���׵���');
