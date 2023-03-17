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
    %%take use of MPNDS to sort
    [FrontNo,MaxFNo] = MPNDS(objs,Population.cons,N,DM);    
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