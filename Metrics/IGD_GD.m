function Score= IGD_GD(PopObj,PF,DM)
% <metric> <min>
%%calculate the IGD and GD for the result
[~,M]=size(PF);
[row,~]=size(PF);
[col,~]=size(PopObj);
Distance_IGD=zeros(row,col);
Distance_GD=zeros(col,row);
%we calculate their distance between real PF and PF for each party, sum
%them,and find the min value for each individual
for i=1:DM
    Distance_IGD = Distance_IGD+pdist2(PF(:,(i-1)*M/DM+1:i*M/DM),PopObj(:,(i-1)*M/DM+1:i*M/DM));
    Distance_GD = Distance_GD+pdist2(PopObj(:,(i-1)*M/DM+1:i*M/DM),PF(:,(i-1)*M/DM+1:i*M/DM));
end  

Distance_IGD = min(Distance_IGD,[],2);
IGD    = mean(Distance_IGD);

Distance_GD = min(Distance_GD,[],2);
GD    = norm(Distance_GD) / length(Distance_GD);

Score=[IGD,GD];
end