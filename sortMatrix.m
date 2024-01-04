function [Q]=sortMatrix(W)
[N,M]=size(W); %计算得出K的二维值
Q=zeros(N,M); %重新定义一个矩阵，将相关性设置成0-1矩阵后放入该矩阵
for i=1:N
    for j=1:M
        if W(i,j)>0
            Q(i,j)=1;
        end
    end
end
end