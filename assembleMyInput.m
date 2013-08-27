function assembleMyInput

dbstop if error;
dbstop if warning;

%rdcFilename1                                   = 'N-H_medium1.m';
%rdcFilename2                                   = 'N-H_medium2.m';
%rdcFilename1                                   = 'N-H_medium.m';
%rdcFilename2                                   = 'C-H_medium.m';
%chRDCs                                        = load('C-H_medium1.m');
combinedResonancesAndProtonCoordinatesFilename = 'combinedResonancesAndProtonCoordinates.txt';
sse_helixFilename                              = 'sse_helix.txt';
%sse_sheetFilename                              = 'sse_sheet.txt';
%nhVectorsFilename                              = 'N-H_vectors.m';
outfilename                                    = 'myinput.m';



%nhRDCs                                         = load(rdcFilename1);
%ccaRDCs                                        = load(rdcFilename2);

[RESNUMS resonanceAA_Name H_CS N_CS protonX protonY protonZ] ...
    = textread(combinedResonancesAndProtonCoordinatesFilename,'%d %s %f %f %f %f %f');


[sseTypes_helix sseIndex_helix sseIndex_helix resName1_helix chainName1_helix resIndex1_helix resName2_helix ...
 chainName2_helix resIndex2_helix someOtherIndex_helix lengthOfSSE_helix] = textread(sse_helixFilename,'%s %d %d %s %s %d %s %s %d %d %d');

%[sheet strandNumberInCurrentSheet sheetId numStrandsInCurrentSheet resName1_sheet chainName1_sheet resIndex1_sheet resName2_sheet chainName2_sheet resIndex2_sheet strandSenseWrtPrevious atomName residueName chainId resSeqNumber atomName residueName chainID resSeqNumber] ...
%    = textread(sse_sheetFilename,'%s %d %s %d %s %s %d %s %s %d %d %s %s %s %d %s %s %s %d');

%nhVectors   = load(nhVectorsFilename);


fprintf(1, 'check out %s\n', outfilename);
fid         = fopen(outfilename, 'w');

for i = 1:length(RESNUMS)
  sseType   = determineSseTypeIsHelixOrNot(RESNUMS(i),...
					   resIndex1_helix,resIndex2_helix);
  
  if (sseType == 'C')
 %   sseType = determineSseTypeIsSheetOrNot(RESNUMS(i),resIndex1_sheet,resIndex2_sheet);
  else
    assert (sseType == 'H');
  end
  
  
  
%  fprintf(fid, '%d\t%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%s\t%s\t%d\t%f\t%f\t%f\n', RESNUMS(i), resonanceAA_Name{i}, nhRDCs(i,2), ccaRDCs(i,2), ...
%		    nhVectors(i,2),nhVectors(i,3),nhVectors(i,4),H_CS(i),N_CS(i),sseType, 'Y',0,protonX(i),protonY(i),protonZ(i));

%  fprintf(fid, '%d\t%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%s\t%s\t%d\t%f\t%f\t%f\n', RESNUMS(i), resonanceAA_Name{i}, nhRDCs(i,2), -999.00, ...
%		    nhVectors(i,2),nhVectors(i,3),nhVectors(i,4),H_CS(i),N_CS(i),sseType, 'Y',0,protonX(i),protonY(i),protonZ(i));

%  fprintf(fid, '%d\t%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%s\t%s\t%d\t%f\t%f\t%f\n', RESNUMS(i), resonanceAA_Name{i}, -999.00, -999.00, ...
%
             %		    nhVectors(i,2),nhVectors(i,3),nhVectors(i,4),H_CS(i),N_CS(i),sseType, 'Y',0,protonX(i),protonY(i),protonZ(i));
	     
	     fprintf(fid, '%d\t%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%s\t%s\t%d\t%f\t%f\t%f\n', RESNUMS(i), resonanceAA_Name{i}, -999.00, -999.00, -999.00,-999.00,-999.00,H_CS(i),N_CS(i),sseType, 'Y',0,protonX(i),protonY(i),protonZ(i));	     

end

fclose(fid);

fprintf(1, 'enter return to exit debug mode.\n');
keyboard
		
%returns 'N' if sseType is determined to be non-helix, 'H' if helix.
function sseType = determineSseTypeIsHelixOrNot(residueNumber, resIndex1, ...
						resIndex2);
		
sseIndex = -1;
assert (length(resIndex1) == length(resIndex2));
for i = 1:length(resIndex1)
  if ((resIndex1(i) <= residueNumber) & (residueNumber <= resIndex2(i)))
   sseIndex = i;
   break;
  end
end

if (sseIndex == -1)
  sseType = 'C';
else
  sseType     = 'H';
%    elseif (strcmp(sseTypeName, 'SHEET') == 1)
%      sseType = 'B';
%    else
%      error('unknown sseType');
%    end
end

%returns 'C' if sseType is determined to be non-sheet, 'B' if sheet.
function sseType = determineSseTypeIsSheetOrNot(residueNumber, resIndex1, ...
						resIndex2);
		

sseType = determineSseTypeIsHelixOrNot(residueNumber, resIndex1, resIndex2);

if (sseType == 'H')
  sseType = 'B';
else
  sseType = 'C';
end

