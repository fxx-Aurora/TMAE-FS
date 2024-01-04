function [pos,neg,bun,Cardinal]=cal_neighbor_low(DC,Neighbor)
[~,D] = size(DC);%决策类别数
[m,~] = size(Neighbor);%样本个数
pos = cell(1,D);
neg = cell(1,D);
bun = cell(1,D);
Upper = cell(1,D);
Cardinal = zeros(1,D);
for d = 1:D
    for i = 1:m
        neighborhood = find((Neighbor(:,i)) == 1);
        if ismember(neighborhood,find(DC(:,d))) == 1
            pos{1,d} = union(pos{1,d},i);
        else
            if isempty(intersect(neighborhood,find(DC(:,d)))) == 1
                neg{1,d} = union(neg{1,d}, i);
            else
                Upper{1,d} = union(Upper{1,d}, i);
            end
        end
    end
    bun{1,d} = setdiff(Upper{1,d},pos{1,d});
    [~,Cardinal(1,d)] = size(pos{1,d});
end
end