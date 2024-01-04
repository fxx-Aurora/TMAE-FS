function [H]=union_matrix(P,L) %对两个矩阵做并集
[N,M]=size(P);
H=zeros(N,M); %定义一个新的矩阵，将两个矩阵中全为1的值进行并集总和到H矩阵中
for i=1:N
    for j=1:M
        if P(i,j)==1||L(i,j)==1
            H(i,j)=1;
        end
    end
end
end