function [Conditional_entropy]=entory(data,rho)  %计算当前数据集下的条件熵
dataD=data(:,end);
dataC=data(:,1:end-1);
[distance,Neighbor,W,~,~]=Neigh(data,rho);
[m,~]=size(dataC);
[DC,DC_diff] = cal_decision(dataD);
[d1,~]=size(DC_diff);
[~, d] = size(DC);

gamma_ds = zeros(2,d);
for j = 1:d
    gamma_ds(1,j) = DC_diff(j);
    A = zeros(m,m);
    A(find(DC(:,j)),find(DC(:,j)))=1;%样本类别为1的基数，行列树为1，其次为0，当j=2
    temp_dis = min(distance,A); %%仅保留同一个类别的距离
    AA = reshape(temp_dis.',1,m*m); 
    AA(find(AA==0)) = [];
    gamma_ds(2,j) = 1/mean(AA); %%相同类训练样本之间距离平均数的倒数
end
%求解当前数据集下的信任函数
[Bel,~,~,~]=DS(Neighbor,dataD,gamma_ds,W,distance);
%信息熵
[N,~]=size(data);  %计算论域data的个数
tmp=0;
for i=1:N
    Neighborhood1=find(Neighbor(:,i)==1);
    [Neighbor_num,~]=size(Neighborhood1);
    target_label=dataD(i); %计算该样本的决策类
    x=0;
    for j=1:Neighbor_num
        if dataD(Neighborhood1(j))==target_label
            x=x+1;  %计算邻域集合的样本决策和该样本决策相同的个数
        end
    end
    for j=1:d1
        if DC_diff(j)==target_label
            B=Bel(i,j);  %样本i对应的信任函数
        end
    end
    tmp=tmp+B*log(x/Neighbor_num);
end
Conditional_entropy=-tmp/N;
end