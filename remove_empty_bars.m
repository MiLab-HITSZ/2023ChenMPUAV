function remove_empty_bars(hBars)
    for iSeries = 1:numel(hBars)
        zData = get(hBars(iSeries),'ZData'); %# Get the z data
        index = logical(kron(zData(2:6:end,2) == 0,ones(6,1))); %# Find empty bars
        zData(index,:) = nan; %# Set the z data for empty bars to nan
        set(hBars(iSeries),'ZData',zData); %# Update the graphics objects
        if zData<1
            set(hBars(iSeries),'FaceColor','#9FB7A0');
        else
            set(hBars(iSeries),'FaceColor','#36907E');
        end
    end
end