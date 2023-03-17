function OptMPNDS(Global)
% <algorithm> <O>
%OptMPNDS
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------
        Population = Global.Initialization();     
        [~,FrontNo,CrowdDis] = EnvironmentalSelection(Population,Global.N,Global.DM);
        %% Optimization        
        while Global.NotTermination(Population)
            MatingPool = TournamentSelection(2,Global.N,FrontNo,-CrowdDis);
            Offspring  = GA(Population(MatingPool));
            [Population,FrontNo,CrowdDis] = EnvironmentalSelection([Population,Offspring],Global.N,Global.DM);
        end  
end



