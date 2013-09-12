function automatedRunModeModel

%normalModeHD: This script runs HD on 1EF1 homology models. 

dbstop if error
dbstop if warning

modeIndex                   = [7, 8, 9, 10, 11, 78];
modelIndex                  = (1 :11 );
modelIndexBidirectional     = (1 :121 );
  
for i = modeIndex;
    
    if(i == 78)
        for j = modelIndexBidirectional;
            runBetter(i,j);
            fprintf ('Mode %d, Model %d \n', i, j);
        end
    else
        for j = modelIndex;
            fprintf ('Mode %d, Model %d \n', i, j);
            runBetter(i,j);
        end
    end
end




