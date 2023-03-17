function [Population,FrontNo,CrowdDis] = MPSELECT(Population,N,DM)
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
    [FrontNo,~] = NDSort(FrontNo_low,Population.cons,N); 
    %% Population for next generation
    Population = Population(FrontNo==1);
  end