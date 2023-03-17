classdef MPMOP10< PROBLEM
% <problem> <MPMOP>


    methods
        %% Initialization
        function obj = MPMOP10()
           obj.Global.M = 9;
            if isempty(obj.Global.D)
                obj.Global.D = 10;
            end
            obj.Global.lower    = zeros(1,obj.Global.D);
            obj.Global.upper    = ones(1,obj.Global.D);
            obj.Global.DM=3;
            obj.Global.encoding = 'real';
        end
        %% Calculate objective values for each party
        function PopObj = CalObj(obj,PopDec)
            M=obj.Global.M;
            PopObj(:,[1:M/3])=MPMOP_Value('MPMOP5', PopDec(), 0);
            PopObj(:,[M/3+1:2*M/3])=MPMOP_Value('MPMOP5', PopDec(), 1);
            PopObj(:,[2*M/3+1:M])=MPMOP_Value('MPMOP5', PopDec(), 1.5);           
        end
        %% Sample reference points on Pareto front
        function P = PF(obj)
            t1=0;
            t2=1;
            t3=1.5;
            [~,true_y1,true_y2,true_y3]=true_PS(obj.Global.D,'MPMOP10',t1,t2,t3);
            P=[true_y1,true_y2,true_y3];
        end
        
        %% Sample reference points on Pareto optimal set
        function P = PS(obj)
            t1=0;
            t2=1;
            t3=1.5;
            P=true_PS(obj.Global.D,'MPMOP10',t1,t2,t3);
        end
        
    end
end