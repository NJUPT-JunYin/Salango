% function LBCA_trashing0426(~)
function [Bn,Mn] = LBCA(n,n0,t0,O0)
t0
% o0=[1,0,1,0,1,0,0,0,0;
%    1,1,0,1,0,0,0,0,0;
%    0,1,0,1,1,1,1,0,0;
%    0,1,0,1,1,0,1,1,0;
%    1,1,0,0,0,1,1,0,1;];%初始图的矩阵表示
% m=9;%例子中用户设备的个数为9
% n0=5;%小基站个数为5
% o0=[0,0,0,1,1,1,0;
%     0,0,0,1,0,1,0;
%     0,0,0,1,1,0,1;
%     1,1,1,0,0,0,0;
%     1,0,1,0,0,0,0;
%     1,1,0,0,0,0,0;
%     0,0,1,0,0,0,0;];
% m=4;
% n0=3;
%修改之前静态的LBCA算法，加上指数加权平均求节点与基站的距离
%每个时间片中重复调用LBCA算法，其实在每个时间片中都是静态的处理
m=n;%其中m=n是节点数目，n0是基站数目
t=t0;
o0=sparse(O0);
N=m+n0;    %初始图的矩阵的行数,第二个参数为1时表示是行数，其实就是基站与节点数目的和
% E=zeros(N);        %最终画图的矩阵，创建全零矩阵，只需要将有连接的位置改成1，最好还是输出E1，简洁
% flag=zeros(m,N);   %记录是否访问过,除去预分配的基站，只记录第二次分配及以后过程
Bn=zeros(1,m);      %下标，表示哪个用户设备只连接了1个小基站，用来预分配的，记录最终被分配过去的小基站的下标
% F=zeros(m,N);      %9行14列的矩阵，用来存储每个用户设备分配过程的基站下标，小于n0的是小基站，大于n0的是 用户设备
Mn=zeros(1,n0);
C1=zeros(m,1);%记录用户设备是否被分配
% E1=zeros(m,n0);

%处理一下flag中可能出现的误差
% k0=1;
% for i=1:m
%     flag(i,k0+n0)=Inf;
%     k0=k0+1;
% end
%% 开始处理所有用户设备和小基站下标
C=sum(o0,1);   %对矩阵o0的每一行求和，得到每个用户设备连接的小基站个数
A2=zeros(m,n0);   %用来存放每个用户设备连接的小基站的下标
for i=1:m
    for j=1:n0
        if o0(j,i)==1
           A2(i,j)=j;
        end
    end      
end

% B2=zeros(n0,m);   %存放小基站连接的用户设备的下标+N
% for i=1:m
%     for j=1:n0
%         if o0(j,i)==1
%            B2(j,i)=i+N;
%         end
%     end
% end
%到这里为止，step1是完全正确的，预分配处理只连接了一个基站的用户设备时完全正确的


for i=1:1:m
    if C(i)==1%处理只与一个基站相连的用户设备
        for j=1:1:n0
            if A2(i,j)~=0
                Bn(i)=j;
                Mn(j)=Mn(j)+1;
                C1(i)=C1(i)+1;
            end
        end
    else 
        continue;
    end
end
for i=1:1:m
        if C1(i)==0&&C(i)>1%没分配且与多个基站相连接
            J=inf(1,n0);%将某个尚未分配的节点所有可连基站的分配情况抄出来
            for j2=1:n0
                if A2(i,j2)~=0
                    J(1,A2(i,j2))=Mn(1,A2(i,j2));%将Mn中所需要的部分拿出来,个数就是C2对应的值
                end
            end
            min_J=min(J);
            [~,col]=find(J(1,:)==min_J);
            kj=length(find(J(1,:)==min_J));
            if kj==1%最小值只有一个
                Bn(i)=col;
                Mn(col)=Mn(col)+1;
                C1(i)=C1(i)+1;
            else%最小值有多个
                a=col(1);
                Mn(a)=Mn(a)+1;%对基站分配数量的增加
                Bn(i)=a;%记录对应节点分配的基站下标  
                C1(i)=1;
            end
        end
end   
% Bn
% flag
% Mn
% A2
% B2
% C1
% C
% o0
end

% for i=1:m
%     if C(i)==1
%         for j=1:n0
%             if o0(j,i)==1
%                 Bn(i)=j;   %记录只与一个基站相连的用户设备所连接的基站的下标
% %                 P(j)=i;   %记录第一次分配时分给基站的用户设备的下标
%                 Mn(j)=Mn(j)+1;
%                 C1(i)=1;
%                 E1(i,j)=E1(i,j)+1; 
%                 for ii=1:N
%                     flag(i,ii)=1;%若是第一次分配就设置flag中对应用户设备那一行全部为1，表示第一次分配
%                 end
%             end 
%         end
%     end
% end

%% 先将只连接一个基站的用户设备的位置点在E中设为1
% for i=1:m
%     if Bn(i)~=0
%         E(i+n0,Bn(i))=1;
%         E(Bn(i),i+n0)=1;
%     end
% end
%以上在最终画图用的邻接矩阵E中都是正确的
%对于以上的操作E，无论是简单还是复杂矩阵，都是正确的



%% 开始算法的step2
% for i=1:m        %根据之前预分配的结果，得到E1矩阵，记录已经预分配了的结果
%     for j=1:n0
%         E1(i,j)=E(i+n0,j);
%     end
% end
% Mn=sum(E1);         %针对基站，不为0代表已经预分配了用户设备，为0代表之前没有预分配
% C1=sum(E1,2);      %针对用户设备，为1代表已经分配了小基站，为0代表没有分配小基站
% C2=C;              %C记录的原始状态下每个用户设备连接的基站的个数
% for i=1:m
%     if C1(i)==0
%         k1=1;
%         for j=1:n0
%             if A2(i,j)~=0 %将用户设备直接连接的基站入队，基站入队是正确的
%                 F(i,k1)=A2(i,j);
%                 flag(i,j)=flag(i,j)+1;
%                 k1=k1+1;
%             end
%         end
%         for j1=1:N
%             if F(i,j1)<N&&F(i,j1)~=0&&Mn(1,F(i,j1))==0&&C1(i,1)==0%用户设备直接分配给基站进行的处理
%                 Mn(1,F(i,j1))=Mn(1,F(i,j1))+1;
%                 C1(i,1)=C1(i,1)+1;
% %                 E(F(i,j1),i+n0)=E(F(i,j1),i+n0)+1;
% %                 E(i+n0,F(i,j1))=E(i+n0,F(i,j1))+1;
%                 E1(i,F(i,j1))=E1(i,F(i,j1))+1;
%                 Bn(1,i)=Bn(1,i)+F(i,j1);
%                 flag(i,F(i,j1)+n0)=1;
%             else
%                 if F(i,j1)<N&&F(i,j1)~=0&&Mn(1,F(i,j1))~=0&&C1(i,1)==0%这里是扫描到的基站已经有分配的用户设备了，将基站分配的用户设备的下标进站
%                     for k2=1:m
%                         if E1(k2,F(i,j1))~=0&&flag(i,k2+n0)==0&&C(i)<=m
%                             [~,bn]=find(Bn==F(i,j1));
%                             cn=length(find(Bn==F(i,j1)));
%                             for k3=1:cn
%                                 F(i,C(i)+1)=bn(k3)+N;
%                                 C(i)=C(i)+1;%0426修改正确，C记录了F中不为0的个数
%                                 flag(i,bn+n0)=1;
%                             end
%                         end
%                     end
%                 else%0426下面这一段if循环，没有进入,看能不能把它改到下一段中的判断中
%                     %因为不会出现重复进队的情况，所以把条件改成判断F中有没有出现这个下标
%                     if F(i,j1)>N&&flag(i,F(i,j1)-m)==0&&C1(i)==0%C1(i,1)==0&&扫描到的是用户设备，就将原始状态中所有连接的基站的下标进站
%                         for j2=1:n0
%                             bm=find(F(i,j2)~=A2(i,j2));
%                             cm=length(find(F(i,j2)==A2(i,j2)));
%                         end
%                         if bm>0&&cm~=0
%                             for k4=1:cm
%                                 F(i,C(i)+1)=A2(F(i,j1)-N,j2);
%                                 C(i)=C(i)+1;
%                                 flag(i,F(i,j1)-N)=1;
%                             end
%                         end                        
%                     else
%                         if F(i,j1)==0||(F(i,j1)>N&&C1(i,1)==0&&flag(i,F(i,j1)-N)~=0)
%                             break;
%                         end
%                     end
%                 end
%             end 
%         end
%     end  
% end
% 
% for i=1:m
%     if C1(i)==0%连接多个基站的用户设备暂时没有分配基站的时候，进入循环处理
%                %首先检查尚未被分配节点的那条进栈路径，分几种情况讨论（说到底就是最小值的个数问题）
%                %第一种是节点所连接的基站中，存在连接节点数目最小的基站，且未重复
%                %第二种是节点所连接基站中，节点数目最小的基站存在重复
%                %第三种是节点能连基站中，所有基站的连接数目是一样的
%         J=zeros(1,n0);%将某个尚未分配的节点所有可连基站的分配情况抄出来
%         for j1=1:n0
%             J(j1)=1000000;%将不可连基站的连接数目设置为一个极端大的值，不影响后续的求最小值
%         end
%         for j2=1:n0
%             if A2(i,j2)~=0
%                 J(1,A2(i,j2))=Mn(1,A2(i,j2));%将M中所需要的部分拿出来,个数就是C2对应的值
%             end
%         end
%         min_J=min(J);
%         [~,col]=find(J(1,:)==min_J);
%         kj=length(find(J(1,:)==min_J));
%         if kj==1%第一种情况，将节点分配给连接数量最小的基站
%             E1(i,col)=1;%对于E和E1的处理就是分配的过程
% %             E(i,col+n0)=1;
% %             E(i+n0,col)=1;
%             Mn(col)=Mn(col)+1;%对基站分配数量的增加
%             Bn(i)=col;%记录对应节点分配的基站下标
%             C1(i)=1;
%         else 
%             if kj>1&&kj<=C2(i)%第二种情况，将节点分配给下标最小的基站
%                 a=col(1);
%                 E1(i,a)=1;%对于E和E1的处理就是分配的过程
% %                 E(i,a+n0)=1;
% %                 E(i+n0,a)=1;
%                 Mn(a)=Mn(a)+1;%对基站分配数量的增加
%                 Bn(i)=a;%记录对应节点分配的基站下标  
%                 C1(i)=1;
%              end
%          end
%     end
% end
%  
% F
% E1
% Bn
% flag
% Mn
% A2
% B2
% C1
% C2
% C
% o0
% %0426目前来看，LBCA的算法是正确的，只是还没有考虑抖动，可以先将变量输出到主函数中，然后在主函数中处理抖动的情况
% end