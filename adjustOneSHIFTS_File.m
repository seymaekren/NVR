function adjustOneSHIFTS_File(templateModeldataFilename, shiftsFilename)

%based on the template file templateModeldataFilename
%and the SHIFTS file, where
%SHIFTS file is defined below, generates new SHIFTS file
%which contains in the order corresponding to the model data file
%the same information as in the original SHIFTS file.

[vectors,types,baseModelFileResidueIDs,sStruct, H_Bond, allDists,I_AllDists] = loaddata(templateModeldataFilename);

outputShiftsFilename = sprintf('%s.adjusted', shiftsFilename);
				   
shiftsData           = load   (shiftsFilename);
if (length(shiftsData) == 0)
  fprintf(1, 'ERROR: shifts file is empty.\n');
  keyboard
  return;
end
writeShiftsFile               (outputShiftsFilename, shiftsData, ...
			       baseModelFileResidueIDs);

function writeShiftsFile (outputShiftsFilename, shiftsData, baseModelFileResidueIDs)

fout = fopen(outputShiftsFilename, 'w');
fprintf(1, 'writing adjusted SHIFTS data to file %s...\n', outputShiftsFilename);

for i = 1:length(baseModelFileResidueIDs)
  outputIndex = find(shiftsData(:,1) == baseModelFileResidueIDs(i));
  if (size(outputIndex)~= 0)
    fprintf(fout, '%d %f %f\n', shiftsData(outputIndex,1), shiftsData(outputIndex,2), shiftsData(outputIndex,3));
  else
    fprintf(1, 'ERROR: could not find the SHIFT corresponding to residue %d\n', baseModelFileResidueIDs(i));
    fprintf(1, 'remove this line from the model data file, and also regenerate all model files from scratch.\n');
    return;
  end
end

fclose(fout);
