function printBinaryDistances(COLIN, NTH, ALLDISTS, inputModeIndex, inputModelIndex)

filename = sprintf('OutputFiles/1EF1/BinaryDistances_%d_model%d.txt',inputModeIndex, inputModelIndex);
fid = fopen(filename,'w');
fprintf(1, 'printing to %s\n', filename);
for residue1Index = 1:length(COLIN)
    for residue2Index =  1:length(COLIN)
        if (ALLDISTS(COLIN(residue1Index),COLIN(residue2Index)) > NTH)
            fprintf(fid, '0 ');
        else
            fprintf(fid, '1 ');
        end
    end
    fprintf(fid, '\n');
end
fclose(fid);
