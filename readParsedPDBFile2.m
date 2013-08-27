function [CA_Coords, parsedPDBFileResidueIDs] = readParsedPDBFile2(coordsFilename);

fin       = fopen(coordsFilename, 'r');

[readAtomName readResidueName readResidueID readX_Coord readY_Coord readZ_Coord] =textread(coordsFilename,'%s %s %d %f %f %f');
%now checking the contents

numCoords = length(readAtomName);
CA_Coords = zeros(numCoords,3);

parsedPDBFileResidueIDs = zeros(numCoords,1);

i = 1; 
	       
while (i <= numCoords)
  assert (strcmp(readAtomName(i),'CA') == 1); 
				    
  CA_Coords(i,1:3) = [readX_Coord(i) readY_Coord(i) ...
		    readZ_Coord(i)];
  parsedPDBFileResidueIDs (i)  = readResidueID(i);
  i = i + 1;
end

fclose(fin);
%fprintf('Total %d CA coords are read\n', numCoords);

