function printOptimizationFiles (CP, SXCP, SSCP, TP, HDE, refineWithRDCs, ...
    RP1, RP2, NOES, ALLDISTS, COLIN, ...
    NTH, useBayesianMatrix)

VERY_LARGE_NUMBER     = 1E9;
combinedScoringMatrix = zeros(size(CP));

filename = sprintf('OutputFiles/combinedScoringMatrix.txt');
fid = fopen(filename,'w');
fprintf(1, 'printing to %s\n',filename);

for rowIndex = 1:size(combinedScoringMatrix,1)
    for columnIndex = 1:size(combinedScoringMatrix,2)
        if (rowIndex == columnIndex)
            checkEntry(CP(rowIndex,rowIndex),   rowIndex, 'CP');
            checkEntry(SXCP(rowIndex,rowIndex), rowIndex, 'SXCP');
            checkEntry(SSCP(rowIndex, rowIndex),rowIndex, 'SSCP');
            checkEntry(TP(rowIndex,rowIndex),   rowIndex, 'TP');
            checkEntry(HDE(rowIndex,rowIndex),  rowIndex, 'HDE');
            if (refineWithRDCs)
                checkEntry(RP1(rowIndex,rowIndex),  rowIndex, 'RP1');
                checkEntry(RP2(rowIndex,rowIndex),  rowIndex, 'RP2');
            end
        end
        
        if ((CP(rowIndex,columnIndex)   == 0) || ...
                (SXCP(rowIndex,columnIndex) == 0) || ...
                (SSCP(rowIndex,columnIndex) == 0) || ...
                (TP(rowIndex,columnIndex)   == 0) || ...
                (HDE(rowIndex,columnIndex)  == 0))
            
            combinedScoringMatrix(rowIndex,columnIndex) = VERY_LARGE_NUMBER;
        
        elseif (refineWithRDCs && ((RP1(rowIndex,columnIndex) == 0) || (RP2(rowIndex,columnIndex)  == 0)))
            
            combinedScoringMatrix(rowIndex,columnIndex) = VERY_LARGE_NUMBER;
        
        else
            combinedScoringMatrix(rowIndex,columnIndex) = -log(CP(rowIndex,columnIndex)) ...
                -log(TP(rowIndex, columnIndex))...
                -log(HDE(rowIndex,columnIndex))...
                -log(SXCP(rowIndex, columnIndex));
                %-log(SSCP(rowIndex,columnIndex));
                
            
            if (refineWithRDCs)
                combinedScoringMatrix(rowIndex,columnIndex) = ...
                    combinedScoringMatrix(rowIndex, columnIndex) ...
                    - log(RP1(rowIndex,columnIndex)) ...
                    - log(RP2(rowIndex,columnIndex));
	
            end
        end
        
        fprintf(fid, '%f ', combinedScoringMatrix(rowIndex, columnIndex));
        
        if (useBayesianMatrix)
            fprintf(fid, '%f ', bayesianMatrix(rowIndex, columnIndex));
        end
    end
    fprintf(fid, '\n');
end

fclose(fid);
    

printNOE_List(NOES);
printBinaryDistances      (COLIN, NTH, ALLDISTS);
%printBinaryDistancesForQAP(COLIN, NTH, ALLDISTS);
%printDistances(ALLDISTS);


