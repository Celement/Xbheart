clear all;
close all;
points=200000;       level=4;    sr=360; 
%读入ECG信号
load('D:\MATLAB\R2012b\work\ecgdata.mat');
R=data(:,2);
X=R(1:200000);

plot(X(1:points));grid on;axis tight;axis([1,points,-2,5]);
title('ECG信号');
 
swa=zeros(4,points);
swd=zeros(4,points);
signal=X(0*200000+1:1*200000);
 
%算小波系数和尺度系数
for i=1:points-3
  swa(1,i+3)=1/4*signal(i+3-2^0*0)+3/4*signal(i+3-2^0*1)+3/4*signal(i+3-2^0*2)+1/4*signal(i+3-2^0*3);
   swd(1,i+3)=-1/4*signal(i+3-2^0*0)-3/4*signal(i+3-2^0*1)+3/4*signal(i+3-2^0*2)+1/4*signal(i+3-2^0*3);
end
j=2;
while j<=level
   for i=1:points-24
     swa(j,i+24)=1/4*swa(j-1,i+24-2^(j-1)*0)+3/4*swa(j-1,i+24-2^(j-1)*1)+3/4*swa(j-1,i+24-2^(j-1)*2)+1/4*swa(j-1,i+24-2^(j-1)*3);
     swd(j,i+24)=-1/4*swa(j-1,i+24-2^(j-1)*0)-3/4*swa(j-1,i+24-2^(j-1)*1)+3/4*swa(j-1,i+24-2^(j-1)*2)+1/4*swa(j-1,i+24-2^(j-1)*3);
   end
   j=j+1;
end
%画出原信号和尺度系数，小波系数
figure;
subplot(level,1,1); plot(X(1:points)); grid on;axis tight;
title('ECG信号及其在j=1,2,3,4尺度下的尺度系数及小波系数');
for i=1:level
    subplot(level+1,2,2*(i)+1);
    plot(swa(i,:)); axis tight;grid on;xlabel('time');
    ylabel(strcat('a   ',num2str(i)));
    subplot(level+1,2,2*(i)+2);
    plot(swd(i,:)); axis tight;grid on;
    ylabel(strcat('d   ',num2str(i)));
end
 
%画出原图及小波系数
figure;
subplot(level,1,1); plot(real(X(1:points)),'b'); grid on;axis tight;
title('ECG信号及其在j=1,2,3,4尺度下的小波系数');
for i=1:level
    subplot(level+1,1,i+1);
    plot(swd(i,:),'b'); axis tight;grid on;
    ylabel(strcat('d   ',num2str(i)));
end
 
%**************************************求正负极大值对*****************************************%
ddw=zeros(size(swd));
pddw=ddw;
nddw=ddw;
%小波系数的大于0的点
posw=swd.*(swd>0);
%斜率大于0
pdw=((posw(:,1:points-1)-posw(:,2:points))<0);
%正极大值点
pddw(:,2:points-1)=((pdw(:,1:points-2)-pdw(:,2:points-1))>0);
%小波系数小于0的点
negw=swd.*(swd<0);
ndw=((negw(:,1:points-1)-negw(:,2:points))>0);
%负极大值点
nddw(:,2:points-1)=((ndw(:,1:points-2)-ndw(:,2:points-1))>0);
%或运算
ddw=pddw|nddw;
ddw(:,1)=1;
ddw(:,points)=1;
%求出极值点的值,其他点置0
wpeak=ddw.*swd;
wpeak(:,1)=wpeak(:,1)+1e-10;
wpeak(:,points)=wpeak(:,points)+1e-10;
 
%画出各尺度下极值点
figure;
for i=1:level
    subplot(level,1,i);
    plot(wpeak(i,:)); axis tight;grid on;
ylabel(strcat('j=   ',num2str(i)));
end
subplot(4,1,1);
title('ECG信号在j=1,2,3,4尺度下的小波系数的模极大值点');
 
interva2=zeros(1,points);
intervaqs=zeros(1,points);
Mj1=wpeak(1,:);
Mj4=wpeak(3,:);
 
%画出尺度3极值点
figure;
plot (Mj4);
title('尺度3下小波系数的模极大值点');
 
posi=Mj4.*(Mj4>0);
%求正极大值的平均
thposi=(max(posi(1:round(points/4)))+max(posi(round(points/4):2*round(points/4)))+max(posi(2*round(points/4):3*round(points/4)))+max(posi(3*round(points/4):4*round(points/4))))/4;
posi=(posi>thposi/3);
nega=Mj4.*(Mj4<0);
%求负极大值的平均
thnega=(min(nega(1:round(points/4)))+min(nega(round(points/4):2*round(points/4)))+min(nega(2*round(points/4):3*round(points/4)))+min(nega(3*round(points/4):4*round(points/4))))/4;
nega=-1*(nega<thnega/4);
%找出非0点
interva=posi+nega;
loca=find(interva);
for i=1:length(loca)-1
    if abs(loca(i)-loca(i+1))<80
       diff(i)=interva(loca(i))-interva(loca(i+1));
    else
       diff(i)=0;
    end
end
%找出极值对
loca2=find(diff==-2);
%负极大值点
interva2(loca(loca2(1:length(loca2))))=interva(loca(loca2(1:length(loca2))));
%正极大值点
interva2(loca(loca2(1:length(loca2))+1))=interva(loca(loca2(1:length(loca2))+1));
intervaqs(1:points-10)=interva2(11:points);
count=zeros(1,1);
count2=zeros(1,1);
count3=zeros(1,1);
mark1=0;
mark2=0;
mark3=0;
i=1;
j=1;
Rnum=0;
%*************************求正负极值对过零，即R波峰值，并检测出QRS波起点及终点*******************%
while i<points
    if interva2(i)==-1
       mark1=i;
       i=i+1;
       while(i<points&&interva2(i)==0)
          i=i+1;
       end
       mark2=i;
%求极大值对的过零点
       mark3= round((abs(Mj4(mark2))*mark1+mark2*abs(Mj4(mark1)))/(abs(Mj4(mark2))+abs(Mj4(mark1))));
%R波极大值点
       R_result(j)=mark3-10;
       count(mark3-10)=1;
%求出QRS波起点
       kqs=mark3-10;
       markq=0;
     while (kqs>1)&&( markq< 3)
         if Mj1(kqs)~=0
            markq=markq+1;
         end
         kqs= kqs -1;
     end
  count2(kqs)=-1;
  
%求出QRS波终点  
  kqs=mark3-10;
  marks=0;
  while (kqs<points)&&( marks<2)
      if Mj1(kqs)~=0
         marks=marks+1;
      end
      kqs= kqs+1;
  end
  count3(kqs)=-1;
  i=i+60;
  j=j+1;
  Rnum=Rnum+1;
 end
i=i+1;
end
%************************删除多检点，补偿漏检点**************************%
num2=1;
while(num2~=0)
   num2=0;
%j=3,过零点
   R=find(count);
%过零点间隔
   R_R=R(2:length(R))-R(1:length(R)-1);
   RRmean=mean(R_R);
%当两个R波间隔小于0.4RRmean时,去掉值小的R波
for i=2:length(R)
    if (R(i)-R(i-1))<=0.4*RRmean
        num2=num2+1;
        if signal(R(i))>signal(R(i-1))
            count(R(i-1))=0;
        else
            count(R(i))=0;
        end
    end
end
end
 
num1=2;
while(num1>0)
   num1=num1-1;
   R=find(count);
   R_R=R(2:length(R))-R(1:length(R)-1);
   RRmean=mean(R_R);
%当发现R波间隔大于1.6RRmean时,减小阈值,在这一段检测R波
for i=2:length(R)
    if (R(i)-R(i-1))>1.6*RRmean
        Mjadjust=wpeak(4,R(i-1)+80:R(i)-80);
        points2=(R(i)-80)-(R(i-1)+80)+1;
%求正极大值点
        adjustposi=Mjadjust.*(Mjadjust>0);
        adjustposi=(adjustposi>thposi/4);
%求负极大值点
        adjustnega=Mjadjust.*(Mjadjust<0);
        adjustnega=-1*(adjustnega<thnega/5);
%或运算
        interva4=adjustposi+adjustnega;
%找出非0点
        loca3=find(interva4);
        diff2=interva4(loca3(1:length(loca3)-1))-interva4(loca3(2:length(loca3)));
%如果有极大值对,找出极大值对
        loca4=find(diff2==-2);
        interva3=zeros(points2,1)';
        for j=1:length(loca4)
           interva3(loca3(loca4(j)))=interva4(loca3(loca4(j)));
           interva3(loca3(loca4(j)+1))=interva4(loca3(loca4(j)+1));
        end
        mark4=0;
        mark5=0;
        mark6=0;
    while j<points2
         if interva3(j)==-1;
            mark4=j;
            j=j+1;
            while(j<points2&&interva3(j)==0)
                 j=j+1;
            end
            mark5=j;
%求过零点
            mark6= round((abs(Mjadjust(mark5))*mark4+mark5*abs(Mjadjust(mark4)))/(abs(Mjadjust(mark5))+abs(Mjadjust(mark4))));
            count(R(i-1)+80+mark6-10)=1;
            j=j+60;
         end
         j=j+1;
     end
    end
 end
end
%画出原图及标出检测结果
figure;
plot(X(0*200000+1:1*200000)),grid on,axis tight,axis([1,points,-2,5]);
title('ECG信号的R波峰值及QRS波波段');
hold on
plot(count,'r');
plot(count2,'k');
plot(count3,'k');
for i=1:Rnum
    if R_result(i)==0;
        break
    end
    plot(R_result(i),X(R_result(i)),'bo','MarkerSize',10,'MarkerEdgeColor','g');
end
hold off
figure;
for i=1:Rnum-1
    if R_result(i)==0;
        break
    end
H(i)=R_result(i);
Q(i)=R_result(i);
S(i)=R_result(i+1);
W(i)=S(i)-Q(i);
end
plot(H,W);
title('原始RR间期');
xlabel('ms');
ylabel('ms');
hold off
figure;

plot(H/1000,W/1000);
title('换算单位后的的RR间期');
xlabel('s');
ylabel('s');
figure;
D=W/1000
L1=medfilt1(D,100000);
L2=D-L1;
plot(H/1000,L2);
title('去除基线漂移的HRV');
xlabel('s');
ylabel('s');
figure;
C=1./L2;
F=60*C;
plot(H/1000,F);
title('最终的HRV');
xlabel('s');
ylabel('心率');   




















