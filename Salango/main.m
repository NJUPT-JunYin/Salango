clear
clc
%% 图像输入
%切换地图过程中需要更改地图位置，基站位置和通信半径
tic
%南京邮电大学地图包输入
store1=imread('C:\Users\24251\Documents\MATLAB\matlab datawork\南京邮电大学地图包\南京邮电大学图层图.png'); %地图图层导入
store2=imread('C:\Users\24251\Documents\MATLAB\matlab datawork\南京邮电大学地图包\p过的图片.png');   %地图导入
% %温哥华地图包输入
% store1=imread('C:\Users\24251\Documents\MATLAB\matlab datawork\温哥华地图包\温哥华图层图.png');%地图图层导入
% store2=imread('C:\Users\24251\Documents\MATLAB\matlab datawork\温哥华地图包\温哥华地图.png');%地图导入
%计算图像长宽 
[ym,xm]=size(store1);                       %读取图像长宽坐标值（倒过来是因为matlab坐标系与常规不太一致）
%x=x/3;                                     %除以三是因为图像缩放比例为33%
XOY=zeros(xm,ym);                           %XOY为路径识别邻接矩阵
%初始化
m=0.9;                                     %行人总数占总节点数比例_仙林地图
% m=0.1;                                       %行人总数占总节点数比例_温哥华地图
% m0=6;                                      %固定节点的数目0501
n=100;                                      %初始节点总个数为80，其中1~24位行人节点，25~80位车辆节点(先将初始节点改为10个)
% n0=20;                                     %小基站的数目设置为20
n0=20;%修改小基站数量
%x=200;                                     %文件库文件总量
rmax=1;                                    %时间轴长度
%time=10;                                  %显示第time秒时的连通图
R=200;                                     %仙林地图小基站通信半径
% R=250;                                       %温哥华地图小基站通信半径

R0=zeros(n0,n);%0425用来记录距离的矩阵，感觉还是由LBCA算法来实现会比较好
R1=zeros(n0,n);%记录上一个时间片中节点的位置
O0=zeros(n0,n);%记录随每个时间片下节点与基站的连接情况，在通信半径内就设对应位置为1
a=0.9;%加权平均的比例系数
% Et1=zeros(n,n0,rmax-1);%0427这个三维矩阵的每一个二维矩阵记录通过LBCA算法得到的分配结果，但是还没有解决抖动的情况，通过抖动的情况对其进行处理，然后修改这个结果矩阵再画图
St=zeros(rmax,n);%0427记录某个节点分配的基站下标，其实解释将B这个结果存到一起%0501修改列值，左右合并
% At=zeros(n,n0,rmax);%每行代表了对应节点在初始时连接的基站下标就是列的数字
% B2t=zeros(n0,n,rmax);%每一行代表了初始情况下基站连接的节点序号加n0
Mt=zeros(rmax,n0);%每一行中的元素表示对应基站连接的节点数目
change_is=zeros(rmax,1);%记录每个时间片中出现节点变动的基站个数
trashing_is=zeros(rmax,1);%和change_is类似的列向量，记录的是发生抖动的基站的个数
Bnt=zeros(rmax,n);%0514记录LBCA算法的结果
Mnt=zeros(rmax,n0);%记录未抖动情况下基站负载情况
% Mnt0=zeros(rmax,n0);%记录常规算法的基站负载情况
% Sn=zeros(rmax,n);%0513调用未考虑抖动的算法，得到未考虑抖动情况下的分配过程
change_isn=zeros(rmax,1);%记录为考虑抖动算法中变更基站的个数

%% ---------------------------图像识别模块----------------------------------
%------------------------------------------------------------------------ 
%XOY中存储着当前节点的移动路径方向信息


% -----------------------读取红色节点模块-------------------------------
% 所要查找的红色点的r, g, b值
r_value=255; 
g_value=0;
b_value=0;
% image 的r, g, b三个分量图像
r=store1(:,:,1);
g=store1(:,:,2);
b=store1(:,:,3);
% 标示出图像image中红色点的位置为1，其它点为0，结果存放在index中
index_r=(r==r_value);
index_g=(g==g_value);
index_b=(b==b_value);
index=index_r&index_g&index_b;
% 最终的红色点位置(x, y)坐标
[yr,xr]=find(index==1);
length=size(xr);
for i=1:1:length(1)
XOY(xr(i),yr(i))=1;       %红色区域路径均置为1
end
%------------------------------------------------------------------


%-----------------------读取蓝色节点模块-------------------------------
% 所要查找的蓝色点的r, g, b值
r_value=0; 
g_value=0;
b_value=255;
% image 的r, g, b三个分量图像
r=store1(:,:,1);
g=store1(:,:,2);
b=store1(:,:,3);
% 标示出图像image中蓝色点的位置为2，其它点为0，结果存放在index中
index_r=(r==r_value);
index_g=(g==g_value);
index_b=(b==b_value);
index=index_r&index_g&index_b;
% 最终的蓝色点位置(x, y)坐标
[yb,xb]=find(index==1);
length=size(xb);
for i=1:1:length(1)
XOY(xb(i),yb(i))=2;       %蓝色区域路径均置为2
end
%------------------------------------------------------------------


% -----------------------读取绿色节点模块-------------------------------
% 所要查找的绿色点的r, g, b值
r_value=0; 
g_value=255;
b_value=0;
% image 的r, g, b三个分量图像
r=store1(:,:,1);
g=store1(:,:,2);
b=store1(:,:,3);
% 标示出图像image中绿色点的位置为3，其它点为0，结果存放在index中
index_r=(r==r_value);
index_g=(g==g_value);
index_b=(b==b_value);
index=index_r&index_g&index_b;
% 最终的绿色点位置(x, y)坐标
[yg,xg]=find(index==1);
length=size(xg);
for i=1:1:length(1)
XOY(xg(i),yg(i))=3;       %绿色区域路径均置为3
end
%------------------------------------------------------------------


% -----------------------读取白色节点模块-------------------------------
% 所要查找的白色点的r, g, b值
r_value=255; 
g_value=255;
b_value=255;
% image 的r, g, b三个分量图像
r=store1(:,:,1);
g=store1(:,:,2);
b=store1(:,:,3);
% 标示出图像image中白色点的位置为1，其它点为0，结果存放在index中
index_r=(r==r_value);
index_g=(g==g_value);
index_b=(b==b_value);
index=index_r&index_g&index_b;
% 最终的白色点位置(x, y)坐标
[yw,xw]=find(index==1);
length=size(xw);
for i=1:1:length(1)
XOY(xw(i),yw(i))=4;       %白色区域路径均置为4
end
%------------------------------------------------------------------


% -----------------------读取黄色节点模块-------------------------------
% 所要查找的黄色点的r, g, b值
r_value=255; 
g_value=255;
b_value=0;
% image 的r, g, b三个分量图像
r=store1(:,:,1);
g=store1(:,:,2);
b=store1(:,:,3);
% 标示出图像image中黄色点的位置为100，其它点为0，结果存放在index中
index_r=(r==r_value);
index_g=(g==g_value);
index_b=(b==b_value);
index=index_r&index_g&index_b;
% 最终的黄色点位置(x, y)坐标
[yy,xy]=find(index==1);
length=size(xy);
for i=1:1:length(1)
XOY(xy(i),yy(i))=-1;       %黄色区域路径均置为-1
end
%------------------------------------------------------------------


% -----------------------读取紫色节点模块-------------------------------
% 所要查找的紫色点的r, g, b值
r_value=255; 
g_value=0;
b_value=255;
% image 的r, g, b三个分量图像
r=store1(:,:,1);
g=store1(:,:,2);
b=store1(:,:,3);
% 标示出图像image中紫色点的位置为5，其它点为0，结果存放在index中
index_r=(r==r_value);
index_g=(g==g_value);
index_b=(b==b_value);
index=index_r&index_g&index_b;
% 最终的紫色点位置(x, y)坐标
[yp,xp]=find(index==1);
length=size(xp);
for i=1:1:length(1)
XOY(xp(i),yp(i))=5;       %紫色区域路径均置为5
end
%------------------------------------------------------------------


%% -----------------------读取转弯一节点模块-------------------------------
% 所要查找的转弯一点的r, g, b值
r_value=1; 
g_value=255;
b_value=255;
% image 的r, g, b三个分量图像
r=store1(:,:,1);
g=store1(:,:,2);
b=store1(:,:,3);
% 标示出图像image中转弯一点的位置，结果存放在index中
index_r=(r==r_value);
index_g=(g==g_value);
index_b=(b==b_value);
index=index_r&index_g&index_b;
% 最终的转弯一点位置(x, y)坐标
[yt1,xt1]=find(index==1);
length=size(xt1,1);
for i=1:1:length
XOY(xt1(i),yt1(i))=1.255;       %转弯一区域路径均置为1.255
end
%------------------------------------------------------------------


%% -----------------------读取转弯二节点模块-------------------------------
% 所要查找的转弯二点的r, g, b值
r_value=2; 
g_value=255;
b_value=255;
% image 的r, g, b三个分量图像
r=store1(:,:,1);
g=store1(:,:,2);
b=store1(:,:,3);
% 标示出图像image中转弯二点的位置，结果存放在index中
index_r=(r==r_value);
index_g=(g==g_value);
index_b=(b==b_value);
index=index_r&index_g&index_b;
% 最终的转弯二点位置(x, y)坐标
[yt2,xt2]=find(index==1);
length=size(xt2,1);
for i=1:1:length
XOY(xt2(i),yt2(i))=2.255;       %转弯二区域路径均置为2.255
end
%------------------------------------------------------------------


%% -----------------------读取转弯三节点模块-------------------------------
% 所要查找的转弯三点的r, g, b值
r_value=3; 
g_value=255;
b_value=255;
% image 的r, g, b三个分量图像
r=store1(:,:,1);
g=store1(:,:,2);
b=store1(:,:,3);
% 标示出图像image中转弯三点的位置，结果存放在index中
index_r=(r==r_value);
index_g=(g==g_value);
index_b=(b==b_value);
index=index_r&index_g&index_b;
% 最终的转弯三点位置(x, y)坐标
[yt3,xt3]=find(index==1);
length=size(xt3,1);
for i=1:1:length
XOY(xt3(i),yt3(i))=3.255;       %转弯三区域路径均置为3.255
end
%------------------------------------------------------------------


%% -----------------------读取转弯四节点模块-------------------------------
% 所要查找的转弯四点的r, g, b值
r_value=4; 
g_value=255;
b_value=255;
% image 的r, g, b三个分量图像
r=store1(:,:,1);
g=store1(:,:,2);
b=store1(:,:,3);
% 标示出图像image中转弯四点的位置，结果存放在index中
index_r=(r==r_value);
index_g=(g==g_value);
index_b=(b==b_value);
index=index_r&index_g&index_b;
% 最终的转弯四点位置(x, y)坐标
[yt4,xt4]=find(index==1);
length=size(xt4,1);
for i=1:1:length
XOY(xt4(i),yt4(i))=4.255;       %转弯四区域路径均置为4.255
end
%------------------------------------------------------------------


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%节点创建模块%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %---------------------------固定节点生成模块-------------------------
%---------------------------小基站节点生成-------------------------------
% S(1).signal='*';%一开始8个节点固定的坐标
% S(1).xd=600;                           
% S(1).yd=200;                            
% S(2).signal='*';
% S(2).xd=250;                           
% S(2).yd=400;
% S(3).signal='*';
% S(3).xd=500;                           
% S(3).yd=600;
% S(4).signal='*';
% S(4).xd=700;                           
% S(4).yd=450;
% S(5).signal='*';
% S(5).xd=400;                           
% S(5).yd=800;                            
% S(6).signal='*';
% S(6).xd=550;                           
% S(6).yd=1000;
% S(7).signal='*';
% S(7).xd=700;                           
% S(7).yd=1200;
% S(8).signal='*';
% S(8).xd=400;                           
% S(8).yd=1300;%一开始固定的8个基站的位置
i1=1;
for x1=400:300:700%仙林地图上基站的位置
    for y1=200:200:1200
        S(i1).xd=x1;
        S(i1).yd=y1;
        S(i1).signal='*';
        i1=i1+1;
    end
end
S(13).xd=600;
S(13).yd=100;
S(13).signal='*';
S(14).xd=250;
S(14).yd=300;
S(14).signal='*';
S(15).xd=200;
S(15).yd=500;
S(15).signal='*';
S(16).xd=100;
S(16).yd=400;
S(16).signal='*';
S(17).xd=500;
S(17).yd=500;
S(17).signal='*';
S(18).xd=500;
S(18).yd=1000;
S(18).signal='*';
S(19).xd=800;
S(19).yd=500;
S(19).signal='*';
S(20).xd=400;
S(20).yd=1350;
S(20).signal='*';%仙林地图上基站的位置坐标

%更换小基站数量时小基站数目的设定,小基站初始设置为基本全覆盖,差不多是10个
% S(1).xd=350;
% S(1).yd=1350;
% S(1).signal='*';
% S(2).xd=690;
% S(2).yd=1250;
% S(2).signal='*';
% S(3).xd=450;
% S(3).yd=1020;
% S(3).signal='*';
% S(4).xd=740;
% S(4).yd=870;
% S(4).signal='*';
% S(5).xd=340;
% S(5).yd=700;
% S(5).signal='*';
% S(6).xd=625;
% S(6).yd=650;
% S(6).signal='*';
% S(7).xd=150;
% S(7).yd=400;
% S(7).signal='*';
% S(8).xd=750;
% S(8).yd=400;
% S(8).signal='*';
% S(9).xd=400;
% S(9).yd=350;
% S(9).signal='*';
% S(10).xd=600;
% S(10).yd=100;
% S(10).signal='*';

%增加小基站的个数为15个
% for x1=400:300:700%仙林地图上基站的位置
%     for y1=200:200:1200
%         S(i1).xd=x1;
%         S(i1).yd=y1;
%         S(i1).signal='*';
%         i1=i1+1;
%     end
% end
% S(13).xd=400;
% S(13).yd=1400;
% S(13).signal='*';
% S(14).xd=250;
% S(14).yd=550;
% S(14).signal='*';
% S(15).xd=150;
% S(15).yd=300;
% S(15).signal='*';

%小基站数量为25
% for x1=400:150:700%仙林地图上基站的位置
%     for y1=200:200:1200
%         S(i1).xd=x1;
%         S(i1).yd=y1;
%         S(i1).signal='*';
%         i1=i1+1;
%     end
% end
% S(19).xd=350;
% S(19).yd=1300;
% S(19).signal='*';
% S(20).xd=250;
% S(20).yd=300;
% S(20).signal='*';
% S(21).xd=200;
% S(21).yd=550;
% S(21).signal='*';
% S(22).xd=100;
% S(22).yd=400;
% S(22).signal='*';
% S(23).xd=500;
% S(23).yd=500;
% S(23).signal='*';
% S(24).xd=300;
% S(24).yd=800;
% S(24).signal='*';
% S(25).xd=800;
% S(25).yd=500;
% S(25).signal='*';


% for x1=100:200:1100%温哥华地图上基站的位置0617
%     for y1=100:200:500
%         S(i1).xd=x1;
%         S(i1).yd=y1;
%         S(i1).signal='*';
%         i1=i1+1;
%     end
% end
% S(19).xd=1200;
% S(19).yd=400;
% S(19).signal='*';
% S(20).xd=1200;
% S(20).yd=600;
% S(20).signal='*';%温哥华地图上基站的位置坐标
for i=1:1:rmax
    for j=1:1:n0
        temp_rnd0=j;
        SC(i,j).xd=S(j).xd;                  %将当前时间轴各位节点位置保存到二维数组中
        SC(i,j).yd=S(j).yd;
        SC(i,j).signal=S(j).signal;     
    end
end

%% ----------------------- 节点记录与移动模块---------------------------------
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%节点创建模块%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ---------------------------节点生成模块-------------------------
    for i=1:1:n
        temp_rnd0=i;
        %---------------------------行人节点生成-------------------------------
        if (temp_rnd0<m*n+1)
            S(i).signal='+';
            origin=fix(size(xp,1)*rand(1));
            while(origin==0)
                origin=fix(size(xp,1)*rand(1));
            end
            S(i).xd=xp(origin);
            S(i).yd=yp(origin);
            S(i).whichcell=0;
        end
        %---------------------------车辆节点生成--------------------------------
        if(temp_rnd0>=m*n+1)%车辆是移动节点                      
            S(i).signal='.';
            S(i).whichcell=0;
            RANDOM_possibility=rand(1);
            if(RANDOM_possibility<=0.25)
                RANDOM_number=1;                              %在红色路径上生成节点
            end
            if(RANDOM_possibility>0.25&&RANDOM_possibility<=0.5)
                RANDOM_number=2;                              %在蓝色路径上生成节点
            end                                                  
            if(RANDOM_possibility>0.5&&RANDOM_possibility<=0.75)
                RANDOM_number=3;                              %在绿色路径上生成节点
            end
            if(RANDOM_possibility>0.75&&RANDOM_possibility<=1.0)
                RANDOM_number=4;                              %在白色路径上生成节点
            end
            switch(RANDOM_number)                                 %记录节点位置
                case 1
                    origin=fix(size(xr,1)*rand(1));
                    while(origin==0)
                        origin=fix(size(xr,1)*rand(1));
                    end
                    S(i).xd=xr(origin);
                    S(i).yd=yr(origin);
                case 2
                    origin=fix(size(xb,1)*rand(1));
                    while(origin==0)
                        origin=fix(size(xb,1)*rand(1));
                    end
                    S(i).xd=xb(origin);
                    S(i).yd=yb(origin);
                case 3
                    origin=fix(size(xg,1)*rand(1));
                    while(origin==0)
                        origin=fix(size(xg,1)*rand(1));
                    end
                    S(i).xd=xg(origin);
                    S(i).yd=yg(origin);
                case 4
                    origin=fix(size(xw,1)*rand(1));
                    while(origin==0)
                        origin=fix(size(xw,1)*rand(1));
                    end
                    S(i).xd=xw(origin);
                    S(i).yd=yw(origin);
                otherwise
                    return
            end
        end     
    end
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%当前路径走向判断模块%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%当前路径走向判断模块%%%%%%%%%%%%%%%%%%%%%%%

for i=1:1:rmax
    for j=1:1:n
        temp_rnd0=j;
       %------------------------人的移动----------------------------------
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if(temp_rnd0<m*n+1)
            RANDOM_possibility=rand(1,1);           %产生随机数
            %0-0.25的概率间命名为1事件，0.25-0.5的概率间命名为2事件，0.5-0.75的概率间命名为3事件，0.75-1的概率间命名为4事件，
             if(RANDOM_possibility<=0.25)
                    RANDOM_number=1;
             end
             if(RANDOM_possibility>0.25&&RANDOM_possibility<=0.5)
                    RANDOM_number=2;
             end
             if(RANDOM_possibility>0.5&&RANDOM_possibility<=0.75)
                    RANDOM_number=3;
             end
             if(RANDOM_possibility>0.75&&RANDOM_possibility<=1.0)
                    RANDOM_number=4;
             end
        %根据上述的随机数产生的事件来决定最后要怎么移动，达成随机移动的目的
        if (temp_rnd0<m*n+1)
            S(j).signal='+';                %以+号表示行人，以下是行人的移动，行人可以出现在地图上的任何地方
           switch(RANDOM_number)
               case 1                       %1情况向上移动一单元
                   S(j).xd=S(j).xd;
                   S(j).yd=S(j).yd+1;
               case 2                       %2情况向下移动一单元
                   S(j).xd=S(j).xd;
                   S(j).yd=S(j).yd-1;
               case 3                       %3情况向左移动一单元
                   S(j).xd=S(j).xd-1;
                   S(j).yd=S(j).yd;               
               case 4                       %4情况向右移动一单元
                   S(j).xd=S(j).xd+1;
                   S(j).yd=S(j).yd;
           end
        end
        %考虑到节点问题，我们将节点固定在地图范围内
        if (S(j).xd>xm||XOY(S(j).xd,S(j).yd)==-1)                       %抵达右界或黄色右界，无法通过
            S(j).xd=S(j).xd-1;
        end
        if (S(j).xd<0||XOY(S(j).xd,S(j).yd)==-1)                        %抵达左界，无法通过
            S(j).xd=S(j).xd+1;
        end
        if (S(j).yd>ym||XOY(S(j).xd,S(j).yd)==-1)                       %抵达上界，无法通过
            S(j).yd=S(j).yd-1;
        end 
        if (S(j).yd<0||XOY(S(j).xd,S(j).yd)==-1)                        %抵达下界，无法通过
            S(j).yd=S(j).yd+1;
        end  
        end
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
      
      %-------------------------车辆节点移动------------------------------
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if(temp_rnd0>=m*n+1)
        %-----------------当前节点所在普通位置方向识别---------------------
        if(XOY(S(j).xd,S(j).yd)==1)         %当前路径为红色，向右情况
            %colourcase='right';
            S(j).xd=S(j).xd+1;
        end
        if(XOY(S(j).xd,S(j).yd)==2)     %当前路径为蓝色，向左情况
            %colourcase='left';
            S(j).xd=S(j).xd-1;
        end
        if(XOY(S(j).xd,S(j).yd)==3)     %当前路径为绿色，向下情况
            %colourcase='down';
            S(j).yd=S(j).yd+1;
        end
        if(XOY(S(j).xd,S(j).yd)==4)     %当前路径为白色，向上情况
            %colourcase='up';
            S(j).yd=S(j).yd-1;
        end
        %-------------------当前节点所在转弯位置方向识别--------------------
        if(XOY(S(j).xd,S(j).yd)==1.255)
            %colourcase='up or right';
            if(rand(1,1)>0.5)
               S(j).yd=S(j).yd-1;  %up
            else
               S(j).xd=S(j).xd+1;  %right
            end
        end
        if(XOY(S(j).xd,S(j).yd)==2.255)
            %colourcase='down or right';
            if(rand(1,1)>0.5)
               S(j).yd=S(j).yd+1;  %down
            else
               S(j).xd=S(j).xd+1;  %right
            end
        end
        if(XOY(S(j).xd,S(j).yd)==3.255)
            %colourcase='down or left';
            if(rand(1,1)>0.5)
               S(j).yd=S(j).yd+1;  %down
            else
               S(j).xd=S(j).xd-1;  %left
            end
        end
        if(XOY(S(j).xd,S(j).yd)==4.255)
            %colourcase='up or left';
            if(rand(1,1)>0.5)
               S(j).yd=S(j).yd-1;  %up
            else
               S(j).xd=S(j).xd-1;  %left
            end
        end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %----------------------时间轴记录部分--------------------------
        A(i,j).xd=S(j).xd;                  %将当前时间轴各位节点位置保存到二维数组中
        A(i,j).yd=S(j).yd;
        A(i,j).signal=S(j).signal; 
        
        for k=1:1:n0%求距离矩阵
            if i==1 
                R0(k,j)=sqrt((A(i,j).xd-SC(i,k).xd).^2+(A(i,j).yd-SC(i,k).yd).^2);%初始计算得出的距离
            elseif i>1
                R1(k,j)=R0(k,j);
                R0(k,j)=a*R1(k,j)+(1-a)*sqrt((A(i,j).xd-SC(i,k).xd).^2+(A(i,j).yd-SC(i,k).yd).^2);%加权计算得出的距离
            end
            if R0(k,j)<=R
                O0(k,j)=1;
%                 O0(j+n0,k,i)=1;
            end
        end
    end
    %画图
    figure(1);                         %创建一个窗口1，记录初始的变化情况
    imshow(store2);
    axis on;
    grid;
    hold on;
 
    for j1=1:1:n
        h=plot(A(i,j1).xd,A(i,j1).yd,A(i,j1).signal);      %循环显示n个点
        hold on;                                        %直接在当前图层显示
        text(A(i,j1).xd,A(i,j1).yd,{j1},'FontName','Euclid','fontsize',8);                  %给所有用户节点加上序号标注
        if A(i,j1).signal=='.'
            set(h,'Color','red','LineWidth',1);         %车辆用红色.标记
        end
        if A(i,j1).signal=='+'
            set(h,'Color','blue','LineWidth',1);        %行人用蓝色+标记
        end
    end
    for k0=1:1:n0
        hold on;
        t=plot(SC(i,k0).xd,SC(i,k0).yd,SC(i,k0).signal);
        hold on;
        text(SC(i,k0).xd,SC(i,k0).yd,{k0},'FontName','Euclid','fontsize',8);               %将小基站的序号在图中表示出来
        set(t,'Color','black','LineWidth',1);          %小基站用黑色*标记
        rectangle('Position',[SC(i,k0).xd-R,SC(i,k0).yd-R,2*R,2*R],'Curvature',[1,1],'EdgeColor',[0.3,0.8,0.9]);%0422显示小基站的通信半径
    end
    legend('\fontname{宋体}行人节点','\fontname{宋体}车辆节点','\fontname{宋体}小基站');   
    % set(gca,'FontName','Euclid','fontsize',8);
    t0=i;
    
    [Bn,Mn] = LBCA(n,n0,t0,O0);%未添加抖动的LBCA算法实践
    Bnt=[Bnt;Bn];
    Mnt=[Mnt;Mn];
    
    Mn0 = normal(n,n0,t0,R0);
    Mnt0(t0,:)= Mn0;%记录常规算法情况下基站负载
%     [B,St,M,Mt] = Atlango(n,n0,t0,O0,St,Mt);
%     Mt=[Mt;M];
%     St=[St;B];
%     Mt(t0,:)=M;%0826
% %     St=cat(1,St,B);
% % %     [B,St,B_M] = Atlango_test(n,n0,t0,O0,St);
% % %     Mt(t0,:)=B_M;
%     for j1=1:1:n%这个循环是用来将B写入到St的
%         St(t0,j1)=B(j1);
%     end
% % %     for j=1:1:n
% % %         A(i,j).whichcell=St(i,j);
% % %     end

end

% figure(2)
% axis on;
% grid;
% hold on;
% xlswrite('LBCA',Mnt',10);xlswrite('normal',Mnt0',10);
% boxplot(Mnt');
% figure(3)
% axis on;
% grid;
% hold on;
% sheet=10;
% xlswrite('LBCA',Mnt',sheet);xlswrite('normal',Mnt0',sheet);
%得到St,也就是说是整个分配后的结果（考虑抖动）；此图正确，先注释
%先看有无抖动算法带来的改变分配基站的节点数目对比图
figure(2);%0926 图可以跑完再画
% axis on;
% grid;
% hold on;
% 
% % Change_is=zeros(rmax,1);
% for i=2:1:rmax
% %     Change_is(i)=change_isn(i)-change_is(i);
%     for j=1:1:n
%         if Bnt(i,j)-Bnt(i-1,j)~=0
%             change_isn(i)=change_isn(i)+1;
%         end
%         if St(i,j)-St(i-1,j)~=0
%             change_is(i)=change_is(i)+1;
%         end
%     end
% %     Change_is(i)=change_isn(i)-change_is(i);
% end
% t=1:1:rmax;
% plot(t,change_is/n,'b',t,change_isn/n,'r');
% xlabel('Round');ylabel('Thrashing Ratio');
% legend('Atlango Algorithm','LBCA Algorithm');
toc
%上方的图是正确的
%画出基站的负载均衡方差,箱线图
% Sn_v=std(Mnt,0,1);
% S_v=std(Mt,0,1);
% yyyyyyy=xlsread('Mnt.xls','B137:U140');
% R1=permute(R0,[1,3,2]);
% xlswrite('R0.xls',R1(:,:,1));
%上面这两句是在求指数加权平均的权重，暂时可以先注释掉

% sheet=1;
% xlswrite('W_Mnt_1',Mnt',sheet);xlswrite('W_Mt_1',Mt',sheet);
% sn = Mnt(:)';
% figure(3)
% axis on;grid;hold on;
% boxplot(Mnt);
% xlabel('基站序号');ylabel('负载方差');
% title('未去抖动基站负载的方差示意图');
% 
% figure(4)
% axis on;grid;
% hold on;
% boxplot(Mt);
% xlabel('基站序号');ylabel('负载方差');
% title('去抖动后基站负载的方差示意图');

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%文件模块
% [tuntu,tuntu_bad,A,SC] = fileshare(A,SC,St,rmax,n,n0,Mnt,Mt);
% figure(5)
% axis on;grid;
% hold on;
% llll=1:1:rmax;
% plot(llll,tuntu,'d--b',llll,tuntu_bad,'d--r');
% xlabel('Round');ylabel('Throughput');
% legend('Publish-Subscribe Scheme','Proactive Subscribe Scheme');
% title('吞吐量随轮次的变化情况');
% 
% 
% 
% %%信任值模块
% for i=1:1:rmax
%     for j=1:1:n
%         A(i,j).group=St(i,j);
%     end
% end
% %%交易模块（transaction）
% %%%交易初始化定义
% rm=50;                                       %用户节点的最大通信半径
% rate=0.9;                                     %交易成功率
% % u=0.6;                                        %直接信任值权重
% jl=0.4;                                       %对节点的奖励程度
% cf=0.6;                                       %对节点的惩罚程度
% zh=0.6;                                       %初始信任值
% % eu=0.1;                                       %恶意节点所占比例
% eu=0.1;%0526
% m=n;
% Bili=zeros(5,rmax);
% kkk=1;
% for u=[0.1,0.2,0.3,0.4,0.5,0.6,0.8]
% % for eu=[0.1,0.2,0.3,0.4,0.5,0.7,0.9]
% % for ii=1:1:5
%     ey=eu.*m;  %恶意节点个数
% %     ey=ii*m;
%     randdata=randperm(m);
%     outrand=randdata(1:ey);
%     e=sort(outrand,'ascend');
% 
%     for time=1:1:(rmax+1)
%         for j=1:1:m                                   
%             A(time,j).plist=[];                       %初始化各节点的交易对象列表
%             A(time,j).plist2=[];
%             A(time,j).partnerlist=zeros(m,1);         %初始化可交易对象列表为全零数组
%             A(time,j).Ttransaction=zeros(m,1);        %初始化各个节点与其他节点的交易记录
%             A(time,j).STtransaction=zeros(m,1);       %初始化各个节点与其他节点的交易成功记录
%             A(time,j).DTtransaction=zeros(m,1);       %初始化各个节点与其他节点的交易失败记录
%             A(time,j).IT=zeros(m,1);                  %初始化各个节点对其他节点的直接信任值
%             A(time,j).RT=zeros(m,1);                  %初始化各个节点对其他节点的推荐信任值
%             A(time,j).T=linspace(zh,zh,m);            %初始化各个节点对其他节点的综合信任值
%             A(time,j).Tj=ones(m,1);                   %初始化各个节点的推荐节点
%             A(time,j).Cd=zeros(m,1);                  %初始化各个节点与其他节点的信誉差
%             A(time,j).Rtd=zeros(m,1);                 %初始化各个节点与其他节点的相对信誉差
%             A(time,j).Cre=linspace(0.5,0.5,m);        %初始化各个节点对其他节点的信誉值
%             A(time,j).rate=rate;                      %初始化各个节点交易的成功率
%         end 
%         for j=e
%             A(time,j).rate=0;                         %初始化恶意节点交易的成功率为0
%         end
%     end
% 
% 
% 
%     for time=1:1:rmax
%     %%%%可交易对象（partner）
%         for j=1:1:m
%             for jj=1:1:m
%                 AN=A(time,j).needlist;
%                 AH=A(time,jj).filelist;
%                 C=intersect(AN,AH);
%                 while size(C)>0                       %判断数组直接的需求与拥有关系
%                     A(time,j).partnerlist(jj)=1;      %更新可交易对象列表（time=1）
%                     A(time,j).plist(end+1)=jj;        %将可交易对象的序号存在结构体中
%                     break
%                 end
%             end
%         end
%         % B1=A(time,1).needlist
%         % B2=A(time,1).havelist
%         % B3=A(time,2).needlist
%         % B4=A(time,2).havelist
%         % C=intersect(B1,B4)
%         % C=intersect(B3,B2)
%         % A(time,1).partnerlist
%         % A(time,2).partnerlist
%         % A(time,1).plist
%         % A(time,2).plist
% 
%         %%%%交易对象筛选（所属基站，范围）
%         for j=1:1:m
%             for j2=A(time,j).plist()
%                 r2=sqrt((A(time,j).xd-A(time,j2).xd)^2+(A(time,j).yd-A(time,j2).yd)^2); %计算两节点直接的距离
%                 if ((r2<=rm) || (A(time,j).group==A(time,j2).group))                    %初步筛选 判断是否在通信范围内或是否为同一组（同一基站范围内）
%                     A(time,j).plist2(end+1)=j2;
%                 end
%             end
%         end
% 
%         %%%%确定交易对象（信任值）
%         for j=1:1:m                           %确定每个用户节点当前的交易对象A(time,j).dx
%             if (size(A(time,j).plist2)~= 0)
%                 j4=A(time,j).plist2(1);
%                 xinrenmax=A(time,j).T(j4);
%                 A(time,j).dx=j4;
%                 for j3=A(time,j).plist2()
%                     xinren=A(time,j).T(j3);
%                     if(xinren>xinrenmax)
%                         xinrenmax=xinren;
%                         A(time,j).dx=j3;
%                     end
%                 end
%             end
%             if (size(A(time,j).plist2)==0)
%                 A(time,j).dx=0;              %当A(time,j).dx=0时，用户没有可交易对象，则不进行交易
%             end    
%         end
%         %%%%确定交易的文件
%         for j=1:1:m
%             if(A(time,j).dx~=0)
%                 j1=A(time,j).dx;
%                 B1=A(time,j).needlist;
%                 B2=A(time,j1).filelist;
%                 C=intersect(B1,B2);
%                 A(time,j).wj=C(1);            %默认交易供求关系重合的第一个文件（暂定）
%             end
%             if(A(time,j).dx==0)
%                 A(time,j).wj=0;
%             end
%         end
%         %%%%交易(更新数据)
%         for j=1:1:m
%             A(time+1,j).Ttransaction=A(time,j).Ttransaction;         %初始化各个节点与其他节点的交易记录
%             A(time+1,j).STtransaction=A(time,j).STtransaction;       %初始化各个节点与其他节点的交易成功记录
%             A(time+1,j).DTtransaction=A(time,j).DTtransaction;       %初始化各个节点与其他节点的交易失败记录
%             A(time+1,j).IT=A(time,j).IT;                             %初始化各个节点对其他节点的直接信任值
%             A(time+1,j).RT=A(time,j).RT;                             %初始化各个节点对其他节点的推荐信任值
%             A(time+1,j).T=A(time,j).T;                               %初始化各个节点对其他节点的综合信任值
%             A(time+1,j).Tj= A(time,j).Tj;                            %初始化各个节点的推荐节点
%             A(time+1,j).Cd=A(time,j).Cd;                             %初始化各个节点与其他节点的信誉差
%             A(time+1,j).Rtd=A(time,j).Rtd;                           %初始化各个节点与其他节点的相对信誉差
%             A(time+1,j).Cre=A(time,j).Cre;                           %初始化各个节点对其他节点的信誉值
%             A(time+1,j).rate=A(time,j).rate;                         %初始化各个节点交易的成功率
%         end
% 
%         for j=1:1:m
%             ra=rand(1);                                        
%             if(A(time,j).dx~=0)
%                 num1=A(time,j).dx;
%                 if(ra<=A(time,num1).rate)                                   %交易成功
%                     A(time+1,j).filelist=[A(time,j).filelist A(time,j).wj];  %更新每个用户节点拥有的文件
%                     B=A(time,j).needlist;
%                     x=find(B==A(time,j).wj);
%                     D=A(time,j).needlist;
%                     D(x)=[];
%                     A(time+1,j).needlist=D;                %更新每个用户节点需求的文件
%                     A(time+1,j).Ttransaction(num1)=A(time,j).Ttransaction(num1)+1;        %更新各个节点与其他节点的交易记录
%                     A(time+1,num1).Ttransaction(j)=A(time,num1).Ttransaction(j)+1;
%                     if(A(time,j).rate~=0)                                                     %若请求交易节点为非恶意节点则给出正确评价
%                         A(time+1,j).STtransaction(num1)=A(time,j).STtransaction(num1)+1;      %更新化各个节点与其他节点的交易成功记录
%                         A(time+1,num1).STtransaction(j)=A(time,num1).STtransaction(j)+1;
%                     end
%                     if(A(time,j).rate==0)                                                      %若请求交易节点为恶意节点则给出失败评价
%                          A(time+1,j).STtransaction(num1)=A(time,j).STtransaction(num1)+1;      %更新化各个节点与其他节点的交易失败记录
%                          A(time+1,num1).DTtransaction(j)=A(time,num1).DTtransaction(j)+1;
%                     end
%                 end
%                 if(ra>A(time,num1).rate)                                    %交易失败
%                     A(time+1,j).filelist=A(time,j).filelist;
%                     A(time+1,j).needlist=A(time,j).needlist;
%                     A(time+1,j).Ttransaction(num1)=A(time,j).Ttransaction(num1)+1;        %更新各个节点与其他节点的交易记录
%                     A(time+1,num1).Ttransaction(j)=A(time,num1).Ttransaction(j)+1;
%                     A(time+1,j).DTtransaction(num1)=A(time,j).DTtransaction(num1)+1;      %更新化各个节点与其他节点的交易失败记录
%                     A(time+1,num1).DTtransaction(j)=A(time,num1).DTtransaction(j)+1;
%                 end
%             end
%             if(A(time,j).dx==0)
%                 A(time+1,j).fiellist=A(time,j).filelist;
%                 A(time+1,j).needlist=A(time+1,j).needlist;
%             end
%         end    
% 
%         %%%%更新信任值
%         %%%%%更新各个节点对其他节点的直接信任值
%         for j=1:1:m
%             num2=A(time,j).dx;
%             if(num2~=0)
%                 S=A(time+1,j).STtransaction(num2);
%                 D=A(time+1,j).DTtransaction(num2);
%                 T=A(time+1,j).Ttransaction(num2);
%                 if(S>D)
%                     A(time+1,j).IT(num2)=(S-D)/T;                 
%                 end
%                 if(S<=D)
%                     A(time+1,j).IT(num2)=0;
%                 end
%             end
%         end
% 
%         %%%%%更新各个节点对其他节点的推荐信任值
%         for j=1:1:m
%             for k=1:1:m
%                 AA=A(time+1,k).Ttransaction;
%                 if(any(AA))
%                     fenzi=0;
%                     fenmu=0;
%                     for kk=1:1:m
%                         AAA=A(time+1,k).Ttransaction(kk);
%                         if(AAA>0)
%                             fenzi=fenzi+A(time,j).Cre(kk)*A(time+1,k).IT(j);
%                             fenmu=fenmu+A(time,j).Cre(kk);
%                         end
%                     end
%                     A(time+1,j).RT(k)=fenzi/fenmu;
%                 end
%             end
%         end
% 
%         %%%%%更新各个节点对其他节点的信誉值
%         for j=1:1:m
%             for k=1:1:m
%                 fenzi1=0;
%                 fenmu1=0;
%                 BB=[];
%                 for l=1:1:m
%                     if(A(time,k).Tj(l)==1)
%                         fenzi1=fenzi1+abs(A(time+1,j).RT(k)-A(time+1,l).IT(k));
%                         fenmu1=fenmu1+1;
%                         BIT=A(time+1,l).IT(k);
%                         BB(end+1)=BIT;
%                     end
%                 end
%                 A(time,j).Cd(k)=fenzi1/fenmu1;
%                 Cd=A(time,j).Cd(k);
%                 Std=std(BB);
%                 A(time+1,j).Rtd(k)=Cd/Std;
%                 Rtd=A(time+1,j).Rtd(k);
%                 if(Rtd<1)
%                     A(time+1,j).Cre(k)=A(time,j).Cre(k)+jl*(1-A(time,j).Cre(k))*(1-Rtd);
%                 end
%                 if(Rtd>=1)
%                     A(time+1,j).Cre(k)=A(time,j).Cre(k)-cf*A(time,j).Cre(k)*(1-1/(Rtd));
%                 end
%             end
%         end
% 
%         %%%%%更新各个节点对其他节点的综合信任值
%         for j=1:1:m                                       
%             for k=1:1:m
%                 A(time+1,j).T(k)=A(time,j).T(k);
%                 if(A(time,j).Ttransaction(k)==1)
%                     A(time+1,j).T(k)=u*A(time+1,j).IT(k)+(1-u)*A(time+1,j).RT(k);
%                 end
%             end
%         end
%     %这里开始保存交易成功率0521
%     %为画图加上的数值保存部分
%     transaction=0;
%     stransaction=0;
%     dtransaction=0;
%     for j=1:1:m
%         transaction=transaction+sum(A(time+1,j).Ttransaction);
%         stransaction=stransaction+sum(A(time+1,j).STtransaction);
%         dtransaction=dtransaction+sum(A(time+1,j).DTtransaction);
%     end
%     bili = stransaction/transaction;
%     Bili(kkk,time+1)=bili;   
%     end
%    kkk=kkk+1; 
%     %在这里加入保存的交易成功率，便于后续画图
%     %以上内容为为了画图而加上的数值保存部分
% end
% suctransaction=Bili(1:7,2:rmax+1);
% xlswrite('success_transaction.xls',suctransaction);
% figure(6)
% axis on;grid;
% hold on;
% t=1:1:rmax;
% plot(suctransaction.');
% legend('eu=0.1','eu=0.2','eu=0.3','eu=0.4','eu=0.5','eu=0.7','eu=0.9');
% xlabel('轮次');ylabel('交易成功率');
% title('恶意节点比例对交易成功率的影响');
% suctransaction=Bili(1:7,2:rmax+1);
% xlswrite('u_success_transaction.xls',suctransaction);
% figure(7)
% axis on;grid;
% hold on;
% t=1:1:rmax;
% plot(suctransaction.');
% legend('u=0.1','u=0.2','u=0.3','u=0.4','u=0.5','u=0.6','u=0.8');
% xlabel('轮次');ylabel('交易成功率');
% title('直接信任值对交易成功率的影响');