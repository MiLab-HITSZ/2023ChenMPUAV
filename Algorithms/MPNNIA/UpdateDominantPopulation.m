function [D,FrontNo_low,FrontNo] = UpdateDominantPopulation(D,N,DM)
% Update the dominant population

%------------------------------- Copyright --------------------------------
% Copyright (c) 2022 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------
    
    objs = D.objs;
    [~,uindex,~] = unique(objs,'rows');
    D = D(uindex);
    N = min(N,numel(D));
    objs = D.objs;
    [Nsize,M] = size(objs);
    FrontNo_low = zeros(Nsize,DM);
    for i =1:DM
        FrontNo_low(:,i) = NDSort(objs(:,(i-1)*(M/DM)+1:i*(M/DM)),D.cons,inf);
    end
    [FrontNo,MaxFNo] = NDSort(FrontNo_low,D.cons,N); 
    Next = FrontNo < MaxFNo;
    Last     = find(FrontNo==MaxFNo);
    CrowdDis = CrowdingDistance(objs,FrontNo);
    [~,Rank] = sort(CrowdDis(Last),'descend');
    Next(Last(Rank(1:N-sum(Next)))) = true;
    CrowdDis = CrowdDis(Next);
    D = D(Next);
    FrontNo=FrontNo(Next);
    FrontNo_low=FrontNo_low(Next,:);

    [~,index] = sortrows([FrontNo' -CrowdDis']);
    D = D(index);
    FrontNo = FrontNo(index);
    FrontNo_low = FrontNo_low(index,:);
end