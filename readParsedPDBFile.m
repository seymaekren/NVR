function [N_Coords, H_Coords, parsedPDBFileResidueIDs] = readParsedPDBFile(coordsFilename);

%the pdb file contains H and N coords for each residue, and for
%each structure (there is more than one structure, for each normal
%mode displacement). Some residues may just
%have one of the NH, H coords, in which case this will be handled 
%by discarding the residue(s) for which there is only one of the H
%and N atoms. This is reported by a message.

addpath('/home/home4/apaydin/Mist/Matlab/General');

fin = fopen(coordsFilename, 'r');

[readAtomName readResidueName readResidueID readX_Coord readY_Coord readZ_Coord] =textread(coordsFilename,'%c %s %d %f %f %f');
%now checking the contents

discard = zeros(length(readAtomName),1); %discard(i) = 1 if i'th
					 %line of the
					 %normalModeCoordsFile is to be
					 %discarded.
i = 1; %i : atomIndex
while (i < length(readAtomName)) 
  if (readResidueID(i) ~= readResidueID(i+1))
    discard(i) = 1;
    fprintf(1, 'discarding the line %c %d %f %f %f\n', ...
	    readAtomName(i,1),readResidueID(i,1), readX_Coord(i,1),readY_Coord(i,1),readZ_Coord(i,1));
  else
    i = i + 1; %skipping the next line, since I am sure that won't
               %be discarded either.
  end
  i = i + 1;
end

numDiscarded = sum(discard);
numCoords = (length(readAtomName)-numDiscarded) / 2;
numCoordsCheck = mod((length(readAtomName)-numDiscarded),2); %should
                                                            %be 0
assert (numCoordsCheck == 0);
N_Coords = zeros(numCoords,3);
H_Coords = zeros(numCoords,3);
parsedPDBFileResidueIDs = zeros(numCoords,1);

atomIndex = 0; %corresponds to the output vectors
i = 1;         %corresponds to the input  vectors, should go (at
               %least)  twice as much as atomIndex since the input
               %vectors store H and N in the same array, whereas
               %the output vectors store them separately, and plus
               %there are some discarded lines.
	       
while (i < length(readAtomName))
  if (discard(i) == 0)
    assert (readAtomName(i) == 'N'); %assuming that the parsed PDB
				     %file is ordered such that N is
				     %followed by H.
    atomIndex = atomIndex + 1;
			     
    N_Coords(atomIndex,1:3) = [readX_Coord(i) readY_Coord(i) ...
		    readZ_Coord(i)];
    i = i + 1;
    assert  (discard(i) == 0);
    assert  (readAtomName(i) == 'H');
    assert  (readResidueID(i) == readResidueID(i-1));
    H_Coords(atomIndex,1:3) = [readX_Coord(i) readY_Coord(i) ...
		    readZ_Coord(i)];
    parsedPDBFileResidueIDs (atomIndex)  = readResidueID(i);
  end
  i = i + 1;
end

fclose(fin);
fprintf('Total %d pairs of N,H atoms are read and %d atoms are discarded\n', atomIndex, numDiscarded);
assert (numDiscarded + atomIndex * 2 == length(readAtomName));
