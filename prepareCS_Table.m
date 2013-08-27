function prepareCS_Table
[lineNumbers, residueNumbers, aaNames, atomNames, atomCodes, chemicalShifts,  otherField1, otherField2] = textread('hsri_bmrb.prot','%d %d %s %s %s %f %f %d');

printInOriginalFormat = 0;
printInGuohuiLinsFormat = 1;

if (printInOriginalFormat)
  fid = fopen('hsri_cs.tab','w');
elseif (printInGuohuiLinsFormat)
  fid = fopen('hsri_scoreInputFile.txt','w');
end
%[residueNumbers, chemicalShifts, otherField1, atomNames, otherField2] = textread('ff2_cs_different_format.txt','%d %f %f %s %d');


%fid = fopen('ff2_cs.tab','w');




%printf(fid, '	       N       CA-1          H         CA\n');
if (printInOriginalFormat)
  fprintf(fid, '	       N       CA-1          H       CB-1         CA         CB\n');
elseif (printInGuohuiLinsFormat)
  fprintf(fid, '	      CA      CB              H       N\n');
end

%fprintf(fid, '	       N       CO-1          H       CB-1       CA-1         CB         CA\n');
maxResidueNumber = max(residueNumbers);
minResidueNumber = min(residueNumbers);
numLines         = length(residueNumbers);

for residueIndex = minResidueNumber:maxResidueNumber+1
  line_N    = findLineIndex(residueNumbers,atomNames, residueIndex, 'N');
  line_CA_1 = findLineIndex(residueNumbers, atomNames, residueIndex-1, ...
			    'CA');
  line_CA   = findLineIndex(residueNumbers, atomNames, residueIndex, ...
			    'CA');
%  line_H    = findLineIndex(residueNumbers, atomNames, residueIndex, ...
%			    'HN');
  line_CB_1 = findLineIndex(residueNumbers, atomNames, ...
			    residueIndex-1,'CB');
  line_CB   = findLineIndex(residueNumbers, atomNames, residueIndex, ...
			    'CB');
%  line_CO_1 = findLineIndex(residueNumbers, atomNames, residueIndex-1, ...
%			    'C');
  line_H    = findLineIndex(residueNumbers, atomNames, residueIndex, ...
			  'H');
  cs_N      = findChemicalShift(chemicalShifts, line_N);
  cs_CA_1   = findChemicalShift(chemicalShifts, line_CA_1);
  cs_CA     = findChemicalShift(chemicalShifts, line_CA);
  cs_H      = findChemicalShift(chemicalShifts, line_H);
  cs_CB_1   = findChemicalShift(chemicalShifts, line_CB_1);
  cs_CB     = findChemicalShift(chemicalShifts, line_CB);
%  cs_CO_1   = findChemicalShift(chemicalShifts, line_CO_1);
  %  fprintf(fid, '%5d    %7.3f    %7.3f    %7.3f    %7.3f\n',residueIndex,cs_N, ...

  if (printInOriginalFormat)
    fprintf(fid, '%5d    %7.3f    %7.3f    %7.3f    %7.3f    %7.3f    %7.3f\n',residueIndex,cs_N, ...
	    cs_CA_1, cs_H, cs_CB_1, cs_CA, cs_CB);
  elseif (printInGuohuiLinsFormat)
    fprintf(fid, '%5d    %7.3f    %7.3f    %7.3f    %7.3f\n',residueIndex,cs_CA, cs_CB, cs_H, cs_N);
  end
    

  
  %  fprintf(fid, '%5d    %7.3f    %7.3f    %7.3f    %7.3f    %7.3f    %7.3f    %7.3f\n',720+residueIndex,cs_N, ...
%	  cs_CO_1, cs_H, cs_CB_1, cs_CA_1, cs_CB, cs_CA);
end
fclose(fid);

function lineNumber    = findLineIndex(residueNumbers,atomNames, ...
				   residueIndex, atomName);
lineNumber = -1;
possibleLineNumbers = find(residueNumbers == residueIndex);
for i = 1:length(possibleLineNumbers)
  if (strcmp(atomNames(possibleLineNumbers(i)), atomName))
    lineNumber = possibleLineNumbers(i);
    break;
  end
end
fprintf(1, 'lineNumber = %d\n',lineNumber);
function cs     = findChemicalShift(chemicalShifts, lineNumber);
if (lineNumber == -1)
  cs = 0.0;
else
  cs = chemicalShifts(lineNumber);
end