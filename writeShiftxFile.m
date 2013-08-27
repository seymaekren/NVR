function writeShiftxFile (outputShiftxFilename, resIndex, resName, ...
			  sstruct, cs1,cs2,cs3,cs4,cs5,cs6, baseModelFileResidueIDs, baseModelFileResidueNames);

fout = fopen(outputShiftxFilename, 'w');
fprintf(1, 'writing adjusted SHIFTX data to file %s...\n', outputShiftxFilename);

[averageCS1, averageCS2, averageCS3, averageCS4, averageCS5, averageCS6] = ...
    computeAverageN_And_H_ChemicalShift(cs1,cs2, cs3, cs4, cs5, cs6);

for baseModelIndex = 1:length(baseModelFileResidueIDs)
  outputIndex = find(resIndex(:,1) == baseModelFileResidueIDs(baseModelIndex));
  
  if (size(outputIndex)~= 0)
    fprintf(fout, '%d %s %s %f %f %f %f %f %f\n', resIndex(outputIndex,1), ...
	    resName{outputIndex,1}, sstruct{baseModelIndex,1}, cs1(outputIndex,1),cs2(outputIndex,1),cs3(outputIndex,1),cs4(outputIndex,1),...
	    cs5(outputIndex,1),cs6(outputIndex,1));
  else
    fprintf(1, 'could not find the SHIFT corresponding to residue %d\n', baseModelFileResidueIDs(baseModelIndex));
    fprintf(1, 'will print average values instead \n');
    fprintf(fout, '%d %s %s %f %f %f %f %f %f\n', baseModelFileResidueIDs(baseModelIndex),...
	    baseModelFileResidueNames{baseModelIndex}, ...
	    sstruct{baseModelIndex,1}, averageCS1,...
	    averageCS2,averageCS3,averageCS4,...
	    averageCS5,averageCS6);
  end
end

fclose(fout);



function [averageCS1, averageCS2, averageCS3, averageCS4, averageCS5, averageCS6] = ...
    computeAverageN_And_H_ChemicalShift(cs1,cs2, cs3, cs4, cs5, cs6);

addpath('~/NVR/trunk');

numElements = length(cs1);

assert(numElements == length(cs2));
assert(numElements == length(cs3));
assert(numElements == length(cs4));
assert(numElements == length(cs5));
assert(numElements == length(cs6));


averageCS1 = 0;
averageCS2 = 0;
averageCS3 = 0;
averageCS4 = 0;
averageCS5 = 0;
averageCS6 = 0;

for i=1:numElements
  averageCS1 = averageCS1 + cs1(i);
  averageCS2 = averageCS2 + cs2(i);
  averageCS3 = averageCS3 + cs3(i);
  averageCS4 = averageCS4 + cs4(i);
  averageCS5 = averageCS5 + cs5(i);
  averageCS6 = averageCS6 + cs6(i);
end


  averageCS1 = averageCS1 / numElements;
  averageCS2 = averageCS2 / numElements;
  averageCS3 = averageCS3 / numElements;
  averageCS4 = averageCS4 / numElements;
  averageCS5 = averageCS5 / numElements;
  averageCS6 = averageCS6 / numElements;



