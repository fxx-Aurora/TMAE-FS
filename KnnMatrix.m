function [dist]=KnnMatrix(L,distance)
distances=sort(distance); %排序后的距离
[N,~]=size(distance);
dist=zeros(N,N);
for j=1:N
    neighbor=find(L(:,j)==1);%第j列的为1的个数就是邻域的个数
    [k,~]=size(neighbor);
    for kk=1:k
        for i=1:N
            if distance(i,j)==distances(kk,j)
                dist(i,j)=1;
            end
        end
    end
end