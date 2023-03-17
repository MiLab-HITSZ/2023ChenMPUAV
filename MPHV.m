function result=MPHV(Population,optimum,DM)
[~,M]=size(Population);
result = 1;
    for i=1:DM
        result=min(result,HV(Population(:,(i-1)*(M/DM)+1:i*(M/DM)),optimum(:,(i-1)*(M/DM)+1:i*(M/DM))));
    end
end