function [r]=rho(data)
rh=[1e-5 1e-4 1e-3 1e-2 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9];
[~,N]=size(rh);
for i=1:N
    [~,Neighbor,~,~,~]=Neigh(data,rh(i));
    [NDER_value(i),Neigh_value(i)]=NDER(data,Neighbor);
end
if all(NDER_value(:) == 1)
    NDER_value(:)=0.5;
end
NDER_value=mapminmax(NDER_value,0.001,0.5);
Neigh_value=mapminmax(Neigh_value,0.001,0.5);
for i=1:N
    NDER_Neigh_value(i)=NDER_value(i)+Neigh_value(i);
end
[max_value,max_index]=max(NDER_Neigh_value);
r=rh(max_index);
max_value1=round(max_value,1)+0.1;
end