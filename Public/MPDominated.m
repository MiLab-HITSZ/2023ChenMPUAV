function Dominated=MPDominated(x,y,flaglist)
Dominated = 0;
for i=1:size(flaglist,1)
    Dominated = Dominated | (all(x(flaglist(i,:))<=y(flaglist(i,:))) && any (x(flaglist(i,:))<y(flaglist(i,:))));
    if Dominated
        return
    end
end
end