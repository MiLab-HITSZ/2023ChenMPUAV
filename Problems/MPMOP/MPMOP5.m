classdef MPMOP5< PROBLEM
% <problem> <MPMOP>

    methods
        %% Initialization
        function obj = MPMOP5()
            obj.Global.M = 6;
            if isempty(obj.Global.D)
                obj.Global.D =10;
            end
            obj.Global.lower    = zeros(1,obj.Global.D);
            obj.Global.upper    = ones(1,obj.Global.D);
            obj.Global.DM=2;
            obj.Global.encoding = 'real';
        end
        %% Calculate objective values for each party
        function PopObj = CalObj(obj,PopDec)
            M=obj.Global.M;
            PopObj(:,[1:M/2])=MPMOP_Value('MPMOP5', PopDec(), 0);
            PopObj(:,[M/2+1:M])=MPMOP_Value('MPMOP5', PopDec(),1.5);
        end
        %% Sample reference points on Pareto front
        function P = PF(obj)
            t1=0;
            t2=1.5;
            [~,true_y1,true_y2]=true_PS(obj.Global.D,'MPMOP5',t1,t2);
            P=[true_y1,true_y2];
        end
        
        %% Sample reference points on Pareto optimal set
        function P = PS(obj)
            t1=0;
            t2=1.5;
            P=true_PS(obj.Global.D,'MPMOP5',t1,t2);
        end
        
    end
end