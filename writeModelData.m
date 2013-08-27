function writeModelData (modelDataFilename, NH_BondVectors, H_Coords, parsedPDBFileResidueIDs, ...
			 types, baseModelFileResidueIDs, sStruct, H_Bond);


addpath('/home/home4/apaydin/Mist/Matlab/General');
addpath('/home/home4/apaydin/Mist/NVR/Routines');

%checks the residue Index to make sure we do not put data that does
%not correspond to any existing data line in the basefile.

fout = fopen(modelDataFilename, 'w');
fprintf(1, 'writing model data to file %s...\n', modelDataFilename);


%assuming that residue indices are *not* sorted in the baseModel file.
%assuming they are sorted in the parsed PDB file.
%assuming that there is no missing residues in the parsedPDB files ...

assert (length(parsedPDBFileResidueIDs) >= length(baseModelFileResidueIDs));
assert (length(types) == length(baseModelFileResidueIDs));
assert (length(types) == length(sStruct));
assert (length(types) == length(H_Bond));

parsedPDBFileDisplacementStartIndex = 1; % a block of pdb
                                         % coordinates runs from
                                         % parsedPDBFileDisplacementStartIndex
                                         % .. parsedPDBFileDisplacementEndIndex-1.
parsedPDBFileDisplacementEndIndex   = 1;
parsedPDBFileFirstResidueIndex      = parsedPDBFileResidueIDs(1);

while (parsedPDBFileDisplacementEndIndex <= length(parsedPDBFileResidueIDs))

  %first find the block
  %find the block of lines which delimit the normal mode analysis
  %corresponding to a particular displacement on the NMA.
  
  parsedPDBFileDisplacementStartIndex = parsedPDBFileDisplacementEndIndex;
  parsedPDBFileDisplacementEndIndex   = parsedPDBFileDisplacementEndIndex ...
      + 1;
  while ((parsedPDBFileDisplacementEndIndex <= length(parsedPDBFileResidueIDs))...
	 &...
	 (parsedPDBFileResidueIDs(parsedPDBFileDisplacementEndIndex) ~= ...
	  parsedPDBFileFirstResidueIndex))...
	parsedPDBFileDisplacementEndIndex = parsedPDBFileDisplacementEndIndex + 1;
  end

  %generating one block of modeldata files. separation done
  %elsewhere for the time being.
  
  %     seek the right parsedPDBFileResidueID_Index here;
  
  for baseModelFileResidueID_Index = 1:length(types) 

    %then go according to the base model file to put a
    %line to the o/p file for each baseModelFileResidueIDs.

    baseFileResidueID = baseModelFileResidueIDs(baseModelFileResidueID_Index,1);
    %seek baseFileResidueID in parsedPDBFileResidueIDs(parsedPDBFileDisplacementStartIndex:parsedPDBFileDisplacementEndIndex-1)
    parsedPDBFileResidueID_Index = ...
	find(parsedPDBFileResidueIDs(parsedPDBFileDisplacementStartIndex:parsedPDBFileDisplacementEndIndex-1) == baseFileResidueID) + parsedPDBFileDisplacementStartIndex - 1; %find finds the relative location; making the location absolute with respect to the parsedPDBFileResidueID_Index vector by the addition
    
    assert (parsedPDBFileResidueID_Index >= parsedPDBFileDisplacementStartIndex);
    
    assert (parsedPDBFileResidueIDs(parsedPDBFileResidueID_Index) == baseFileResidueID);
    
%    keyboard
    
    fprintf(fout, '%d\t%s\t%f\t%f\t%f\t%s\t%s\t%f\t%f\t%f\n', ...
	    baseModelFileResidueIDs(baseModelFileResidueID_Index,1), ...
	    types{baseModelFileResidueID_Index,1},...
	    NH_BondVectors(parsedPDBFileResidueID_Index,1), ...
	    NH_BondVectors(parsedPDBFileResidueID_Index,2), ...
	    NH_BondVectors(parsedPDBFileResidueID_Index,3), ...
	    sStruct{baseModelFileResidueID_Index,1}, ...
	    H_Bond{baseModelFileResidueID_Index,1}, ...
	    H_Coords(parsedPDBFileResidueID_Index,1), ...
	    H_Coords(parsedPDBFileResidueID_Index,2), ...
	    H_Coords(parsedPDBFileResidueID_Index,3));
  end
end

fclose(fout);
