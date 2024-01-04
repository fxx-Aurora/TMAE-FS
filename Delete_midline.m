function [Midline,T]=Delete_midline(R,Distance)

[N,~]=size(Distance);
for i=1:N
    %将R矩阵中对应的邻域样本放在Neighbor中
    Neighbor=[];
    for j=1:N
        if R(j,i)==1 && j~=i
           Neighbor=[Neighbor;j];
        end
    end
    %求解每个样本对应的几何矩阵Midline
    [M,~]=size(Neighbor); %得出有多少个邻域样本
    Midline=zeros(M,M); %设置空矩阵
    Midline(1,1)=1;
    for p=2:M  %中垂线的另一端为邻域样本中的第一个样本
        if Distance(Neighbor(p,1),i)>Distance(Neighbor(p,1),Neighbor(1,1))
            Midline(p,1)=2;
        else
            Midline(p,1)=1;
        end
    end
    for k=2:M
        for p=1:M
            if p~=k
                if Distance(Neighbor(p,1),i)>Distance(Neighbor(p,1),Neighbor(k,1))
                    Midline(p,k)=Midline(p,k-1)+1;
                else
                    Midline(p,k)=Midline(p,k-1);
                end
            else
                Midline(p,k)=Midline(p,k-1);
            end
        end
    end
    %取出Midline中的最大值
    %Mid=Midline(:,M);
    max=0;
    for j=1:M
        if max<Midline(j,M)
            max=Midline(j,M);
        end
    end
    
    %将最大值所对应的邻域从邻域集合中删除，即R=0
    if max~=1
        for j=1:M
            if Midline(j,M)==max
                x=Neighbor(j,1);
                R(x,i)=0;
            end
        end
    end
end
T=R;
end