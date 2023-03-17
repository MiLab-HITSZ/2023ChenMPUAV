classdef MPMOP1< PROBLEM
% <problem> <MPMOP>
    methods
        %% Initialization
        function obj = MPMOP1()
            obj.Global.M = 4;
            if isempty(obj.Global.D)
                obj.Global.D =10;
            end
            obj.Global.lower    = zeros(1,obj.Global.D);
            obj.Global.upper    = ones(1,obj.Global.D);
            obj.Global.lower (:,1)=1;
            obj.Global.upper(:,1)=4;
            obj.Global.DM=2;
            obj.Global.encoding = 'real';
        end
        %% Calculate objective values for each party
        function PopObj = CalObj(obj,PopDec)         
            M=obj.Global.M;
            PopObj(:,[1:M/2])=MPMOP_Value('MPMOP1', PopDec(), 1);
            PopObj(:,[M/2+1:M])=MPMOP_Value('MPMOP1', PopDec(), 2);
        end
        %% Sample reference points on Pareto front
        function P = PF(obj)
            t1=1;
            t2=2;
            [~,true_y1,true_y2]=true_PS(obj.Global.D,'MPMOP1',t1,t2);
            P=[true_y1,true_y2];
        end
        
        %% Sample reference points on Pareto optimal set
        function P = PS(obj)
            t1=1;
            t2=2;
            P=true_PS(obj.Global.D,'MPMOP1',t1,t2);
        end
        
    end
end