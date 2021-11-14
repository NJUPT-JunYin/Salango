function pic_change(~)
sheet=1;
A=xlsread('normal.xls',sheet);
% B=xlsread('LBCA.xls',sheet);
C=xlsread('atlango.xls',sheet);
figure(1)
position_1 = 0.8:1:19.8; 
boxplot(A,'colors','b','positions',position_1,'width',0.2,'symbol',''); 
hold on
grid on
position_2 = 1.1:1:20.1; 
% boxplot(B,'colors','r','positions',position_2,'width',0.2,'symbol',''); 
boxplot(C,'colors','r','positions',position_2,'width',0.2,'symbol',''); 
set(gca,'XLim',[0 21]);                          %限定横轴显示长度
set(gca,'YLim',[0 35]);                         %限定纵轴显示长度
xlabel('用户设备数量');ylabel('小基站负载');
hLegend = legend(findall(gca,'Tag','Box'), {'Group A','Group C'});
box_vars = findall(gca,'Tag','Box');
hLegend = legend(box_vars([8,6]), {'常规算法','atlango算法'});
legend('Location','northwest')
% 
% sheet=6;
% At=xlsread('normal_t.xls',sheet);
% Bt=xlsread('LBCA_t.xls',sheet);
% figure(3)
% A1=At(1,:);B1=Bt(1,:);
% n=1:1:20;
% plot(n,A1,'b--','LineWidth',2);
% hold on
% plot(n,B1,'r-','LineWidth',1.5);
% xlabel('小基站');ylabel('小基站负载');
% legend('常规算法','LBCA算法');
% hold on 
% grid on

%0819换一种方式表示基站坐标

end