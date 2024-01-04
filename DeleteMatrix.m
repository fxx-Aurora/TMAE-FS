%根据平均距离，对所选邻居进行删除
function [R]=DeleteMatrix(H,distance)
[~,B]=size(H);
R=H;
for i=1:B
    dist=0;
    neighbor=find(R(:,i)==1);
    [a,~]=size(neighbor);
    for j=1:a
        dist=dist+distance(neighbor(j),i);
    end
    ave=dist/(a-1); %因为邻域集合中包含了自己
    for j=1:a
        if distance(neighbor(j),i)>ave
            R(neighbor(j),i)=0;
        end
    end
end
end