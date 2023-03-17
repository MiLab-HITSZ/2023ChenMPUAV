function [Population,FrontNo,CrowdDis] = EnvironmentalSelection(Population,N,DM)
% The environmental selection of OptMPNDS

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------
    %% Non-dominated sorting
    objs=Population.objs;
    [Nin,M]= size(objs);
    FrontNo_low = zeros(Nin,DM);
    for i =1:DM
        FrontNo_low(:,i) = NDSort(objs(:,(i-1)*(M/DM)+1:i*(M/DM)),Population.cons,Nin);
    end
    [FrontNo,MaxFNo] = NDSort(FrontNo_low,Population.cons,N); 
    Next = FrontNo < MaxFNo;
    CrowdDis = CrowdingDistance(objs,FrontNo);

    Last     = find(FrontNo==MaxFNo);
    [~,Rank] = sort(CrowdDis(Last),'descend');
    Next(Last(Rank(1:N-sum(Next)))) = true;

    %% Population for next generation
    Population = Population(Next);
    FrontNo    = FrontNo(Next);
    CrowdDis   = CrowdDis(Next);
  end