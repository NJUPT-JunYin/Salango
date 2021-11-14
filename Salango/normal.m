%function normal(~)
function Mn0 = normal(n,n0,t0,R0)
t0
%小论文中用于常规算法的对比
% 将用户节点分配给距离最近的小基站，同样是用Mn0记录每个小基站分配的节点数目
Mn0=zeros(1,n0);
% C1=zeros(m,1);
% 取出需要的部分
[~,index]=min(R0,[],1);
Bn0=index';
k=length(index);
for i=1:1:n
    Mn0(index(1,i))=Mn0(index(1,i))+1;
end
end