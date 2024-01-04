function [distance,Neighbor,W,H,R]=Neigh(data,rho)
distance = pdist2(data(:,1:end-1),data(:,1:end-1),'minkowski',2);
FlattenedData = distance(:)'; % 
MappedFlattened = mapminmax(FlattenedData, 0, 1); % 
distance2 = reshape(MappedFlattened, size(distance)); % 
distance = distance2;
%使用lasso形成的邻域矩阵
%rho = 0.01; %系数参数
[W] = leave1(data,rho); 
[Q]=sortMatrix(W); %lasso形成的相关性矩阵中如果该相关性不为0，则视为选中该样本

[~,B]=size(Q);
L=zeros(B,B);
for i=1:B
    for j=1:B
        if i==j
            L(i,j)=1;
        elseif i<j
            L(i,j)=Q(i,j);
        else
            L(i,j)=Q(i-1,j);
        end
    end
end

P=KnnMatrix(L,distance); %knn形成的邻域矩阵
H=union_matrix(P,L); %将两个邻域矩阵做并集
[R]=DeleteMatrix(H,distance); %求出距离的平均值，如果距离大于平均值则将其从邻域集合中删除
[~,Neighbor]=Delete_midline(R,distance);  %根据中垂线得出的距离矩阵，删除具有最大值的邻域，形成的Neighbor矩阵为最终的邻域矩阵，其中邻域包含自己
end