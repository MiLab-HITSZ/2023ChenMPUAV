classdef GLOBAL < handle
%GLOBAL - The class of experimental setting.
%
%   This is the class of experimental setting. An object of GLOBAL class
%   stores the algorithm to be executed, the problem to be solved, and all
%   the parameter settings like population size, number of objectives, and
%   maximum number of evaluations. This class also has several methods
%   which can be invoked by algorithms.
%
% GLOBAL properties:
%   N               <public>	population size
%   M               <read-only>	number of objectives
%   D               <read-only>	number of variables
%   lower           <read-only>	lower bound of each decision variable
%   upper           <read-only>	upper bound of each decision variable
%   algorithm       <read-only>	algorithm function
%   problem         <read-only>	problem function
%   encoding        <read-only> encoding of the problem
%   evaluated       <read-only>	number of evaluated individuals
%   evaluation      <read-only>	maximum number of evaluations
%   gen             <read-only>	current generation
%   maxgen          <read-only>	maximum generation
%   run             <read-only>	run number
%   runtime         <read-only>	runtime
%   save            <read-only> number of saved populations
%   result          <read-only>	set of saved populations
%   PF              <read-only>	true Pareto front
%   parameter       <read-only>	parameters of functions specified by users
%   outputFcn    	<read-only>	function invoked after each generation
%
% GLOBAL methods:
%   GLOBAL          <public>    the constructor, all the properties will be
%                               set when the object is creating
%   Start           <public>    start running the algorithm
%   Initialization  <public>    randomly generate an initial population
%   NotTermination  <public>	terminate the algorithm if the number of
%                               evaluations has exceeded
%   ParameterSet    <public>    obtain the parameter settings from user
%   GetObj          <static>    get the current GLOBAL object

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties
        N          = 100;               % Population size
    end
    properties(SetAccess = ?PROBLEM)
        data;
        M;                              % Number of objectives
        D=10;                              % Number of decision variables
        lower;                          % Lower bound of each decision variable
        upper;                          % Upper bound of each decision variable
        encoding   = 'real';            % Encoding of the problem
        evaluation = 10000;             % Maximum number of evaluations
        DM=1;                           %number of decision makers
        funname={};
        fn;
    end
    properties(SetAccess = ?INDIVIDUAL)
        evaluated  = 0;                 % Number of evaluated individuals
    end
    properties(SetAccess = private)
        algorithm  = @OptAll;       	% Algorithm function
        problem    = @MPMOP1;            % Problem function
        gen;                            % Current generation
        maxgen;                         % Maximum generation
        run        = 1;                 % Run number
        runtime    = 0;                 % Runtime
        save       = 100;             	% Number of saved populations
        result     = {};                % Set of saved populations
        score={}   ;                     %store score
        PF;                             % True Pareto front
        PS;                             %True Pareto optimal set
        parameter  = struct();      	% Parameters of functions specified by users
        outputFcn  = @GLOBAL.Output;  	% Function invoked after each generation
    
    end
    methods
        %% Constructor
        function obj = GLOBAL(varargin)
        %   Example:
        %       GLOBAL('-algorithm',@NSGAII,'-problem',@DTLZ2,'-N',100,...
        %              '-M',2,'-D',10)
            obj.GetObj(obj);
            % Initialize the parameters which can be specified by users
            propertyStr = {'N','M','D','algorithm','problem','evaluation','run','save','outputFcn','rate','data'};
            if nargin > 0
                IsString = find(cellfun(@ischar,varargin(1:end-1))&~cellfun(@isempty,varargin(2:end)));
                [~,Loc]  = ismember(varargin(IsString),cellfun(@(S)['-',S],propertyStr,'UniformOutput',false));
                for i = find(Loc)
                    obj.(varargin{IsString(i)}(2:end)) = varargin{IsString(i)+1};
                end
            end
            % Instantiate a problem object
            if ~isempty(obj.data)
                obj.problem = obj.problem(obj.data);
            else
                obj.problem = obj.problem();
            end
            % Add the folders of the algorithm and problem to the top of
            % the search path
            addpath(fileparts(which(class(obj.problem))));
            addpath(fileparts(which(func2str(obj.algorithm))));
        end
        %% Start running the algorithm
        function Start(obj)
            if obj.evaluated <= 0
                obj.PF = obj.problem.PF();
                obj.PS=obj.problem.PS();
                try
                    tic;
                    obj.algorithm(obj);
                catch err
                    if strcmp(err.identifier,'GLOBAL:Termination')
                        return;
                    else
                        rethrow(err);
                    end
                end
                obj.evaluated = max(obj.evaluated,obj.evaluation);
                if isempty(obj.result)
                    obj.result = {obj.evaluated,INDIVIDUAL()};
                end
            	obj.outputFcn(obj);
            end
        end
        %% Randomly generate an initial population
        function Population = Initialization(obj,N)
        
            if nargin < 2
                N = obj.N;
            end
            Population = INDIVIDUAL(obj.problem.Init(N));
        end
        %% Terminate the algorithm if the number of evaluations has exceeded
        function notermination = NotTermination(obj,Population,Archive)      
            obj.runtime = obj.runtime + toc;
            % Save the population
            if obj.save<=0; num=10; else num=obj.save; end
            index = max(1,min(min(num,size(obj.result,1)+1),ceil(num*obj.evaluated/obj.evaluation)));        
            notermination = obj.evaluated < obj.evaluation;
            if  ~notermination
                if exist('Archive','var')
                    [~,loc,~]=unique(Archive.objs,'rows');
                    Population=Archive(loc);
                end
                obj.result(index,:) = {obj.evaluated,Population};
                Feasible     = find(all(obj.result{end}.cons<=0,2));                             
                objs=obj.result{end}(Feasible).objs;
                NonDominated=true(1,numel(Population));
                %find the solutions that cannot been domianted by other individuals for each party
                for i=1:obj.DM
                    NonDominated =NonDominated & (NDSort(objs(:,(i-1)*obj.M/obj.DM+1:i*obj.M/obj.DM),1) == 1);
                end                
                Population   = obj.result{end}(Feasible(NonDominated));
                Metrics = {@IGD_GD};
                Score     = cellfun(@(S)GLOBAL.Metric(S,Population,obj.PF,obj.DM),Metrics,'UniformOutput',false);  
                %Save the IGD, GD and Solutions Number
                obj.score(index,:)  = {obj.evaluated,Score{1}(1),Score{1}(2),size(Population,2)};                
            end
            assert(notermination,'GLOBAL:Termination','Algorithm has terminated');
            tic;
        end

               
        
        %% Obtain the parameter settings from user
        function varargout = ParameterSet(obj,varargin)
            CallStack = dbstack();
            caller    = CallStack(2).file;
            caller    = caller(1:end-2);
            varargout = varargin;
            if isfield(obj.parameter,caller)
                specified = cellfun(@(S)~isempty(S),obj.parameter.(caller));
                varargout(specified) = obj.parameter.(caller)(specified);
            end
        end
        %% Variable constraint
        function set.N(obj,value)
            obj.Validation(value,'int','size of population ''-N''',1);
            obj.N = value;
        end
        function set.M(obj,value)
            obj.Validation(value,'int','number of objectives ''-M''',2);
            obj.M = value;
        end
        function set.D(obj,value)
            obj.Validation(value,'int','number of variables ''-D''',1);
            obj.D = value;
        end
        function set.algorithm(obj,value)
            if iscell(value) 
                obj.Validation(value{1},'function','algorithm ''-algorithm''');
                obj.algorithm = value{1};
                obj.parameter.(func2str(value{1})) = value(2:end);
            else
                obj.Validation(value,'function','algorithm ''-algorithm''');
                obj.algorithm = value;
            end
        end
        function set.problem(obj,value)
            if iscell(value)
                obj.Validation(value{1},'function','test problem ''-problem''');
                obj.problem = value{1};
                obj.parameter.(func2str(value{1})) = value(2:end);
            elseif ~isa(value,'PROBLEM')
                obj.Validation(value,'function','test problem ''-problem''');
                obj.problem = value;
            else
                obj.problem = value;
            end
        end
        function set.evaluation(obj,value)
            obj.Validation(value,'int','number of evaluations ''-evaluation''',1);
            obj.evaluation = value;
        end
        function set.run(obj,value)
            obj.Validation(value,'int','run number ''-run''',1);
            obj.run = value;
        end
        function set.save(obj,value)
            obj.Validation(value,'int','number of saved populations ''-save''',0);
            obj.save = value;
        end
        %% Variable dependence
        function value = get.gen(obj)
            value = ceil(obj.evaluated/obj.N);
        end
        function value = get.maxgen(obj)
            value = ceil(obj.evaluation/obj.N);
        end
    end
    methods(Static)
        %% Get the current GLOBAL object
        function obj = GetObj(obj)       
            persistent Global;
            if nargin > 0
                Global = obj;
            else
                obj = Global;
            end
        end        
    end
    

    % The following functions cannot be invoked by users
    methods(Access = private)
        %% Check the validity of the specific variable
        function Validation(obj,value,Type,str,varargin)
            switch Type
                case 'function'
                    assert(isa(value,'function_handle'),'INPUT ERROR: the %s must be a function handle',str);
                    assert(~isempty(which(func2str(value))),'INPUT ERROR: the function <%s> does not exist',func2str(value));
                case 'int'
                    assert(isa(value,'double') && isreal(value) && isscalar(value) && value==fix(value),'INPUT ERROR: the %s must be an integer scalar',str);
                    if ~isempty(varargin); assert(value>=varargin{1},'INPUT ERROR: the %s must be not less than %d',str,varargin{1}); end
                    if length(varargin) > 1; assert(value<=varargin{2},'INPUT ERROR: the %s must be not more than %d',str,varargin{2}); end
                    if length(varargin) > 2; assert(mod(value,varargin{3})==0,'INPUT ERROR: the %s must be a multiple of %d',str,varargin{3}); end
            end
        end
    end
    methods(Access = private, Static)       
        function cb_metric(hObject,eventdata,obj,metric)
            metricName   = func2str(metric);
            MetricValues = get(gcbf,'UserData');
            % Calculate the specified metric value of each population
            if ~isfield(MetricValues,func2str(metric)) 
                tempText = text('Units','normalized','Position',[.4 .5 0],'String','Please wait ... ...'); drawnow();
                MetricValues.(metricName)(:,1) = obj.result(:,1);
                MetricValues.(metricName)(:,2) = cellfun(@(S)GLOBAL.Metric(metric,S,obj.PF),obj.result(:,2),'UniformOutput',false);
                set(gcbf,'UserData',MetricValues);
                delete(tempText);
            end
            % Display the results
            cla; Draw(cell2mat(MetricValues.(metricName)),'-k.','LineWidth',1.5,'MarkerSize',10);
            xlabel('Number of Evaluations');
            ylabel(metricName);
            GLOBAL.cb_menu(hObject);
        end
        function cb_menu(hObject)
            % Switch the selected menu
            set(get(get(hObject,'Parent'),'Children'),'Checked','off');
            set(hObject,'Checked','on');
        end
        
        function value = Metric(metric,Population,PF,DM)
            try
                value = metric(Population.objs,PF,DM);
            catch
                value = [inf inf];
            end
        end
        
    end
end