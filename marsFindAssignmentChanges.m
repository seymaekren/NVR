function marsFindAssignmentChanges(includeUnidirectionalMode,includeBidirectionalMode)

addpath('/home/home4/apaydin/Mist/Matlab/Other');
addpath('/home/home4/apaydin/Mist/Matlab/General');
addpath('/home/home4/apaydin/Mist/NVR/Routines');

if (nargin == 0)
  includeUnidirectionalMode = 1;
  includeBidirectionalMode  = 1;
elseif (nargin == 1)
  assert (includeUnidirectionalMode == 1);
  includeBidirectionalMode  = 0;
end


numModelsPerMode       = 11;
highestNormalModeIndex = 11;

numBidirectionalModels = 121;

dbstop if error;

%%
if (includeBidirectionalMode)

    for model1Index = 1:numBidirectionalModels-1

      model2Index                   = model1Index+1;

      assignments1Filename          = sprintf('Bimodal/Results7.8_%d/confidentlyAssignedResidues.txt',model1Index);

      assignments2Filename          = sprintf('Bimodal/Results7.8_%d/confidentlyAssignedResidues.txt',model2Index);

      if (~checkIfFileExists(assignments1Filename))
        fprintf(1, 'couldnt read %s\n',assignments1Filename);
        keyboard
        continue;
      end

      if (~checkIfFileExists(assignments2Filename))
        fprintf(1, 'couldnt read %s\n',assignments2Filename);
        keyboard
        continue;
      end

      assignments1         = load (assignments1Filename);

      assignments2         = load (assignments2Filename);


      assignments2         = compareAssignments1And2(assignments1,assignments2,model1Index,model2Index);

      if (~isempty(assignments2))
        compareAssignments1And2(assignments2,assignments1,model2Index,model1Index);
      end

      fprintf(1, '-----------------------------------------------------------\n');
      fprintf(1, '-----------------------------------------------------------\n');
    end
end
  
%%
if (includeUnidirectionalMode)
    for modeIndex = 7:highestNormalModeIndex
      for modelIndex = 1:numModelsPerMode-1

        assignments1Filename          = sprintf('Unidirectional/Results%d_%d/confidentlyAssignedResidues.txt',modeIndex,modelIndex);
        assignments2Filename          = sprintf('Unidirectional/Results%d_%d/confidentlyAssignedResidues.txt',modeIndex,modelIndex+1);

          if (~checkIfFileExists(assignments1Filename))
            fprintf(1, 'couldnt read %s\n',assignments1Filename);
            keyboard
            continue;
          end

          if (~checkIfFileExists(assignments2Filename))
            fprintf(1, 'couldnt read %s\n',assignments2Filename);
            keyboard
            continue;
          end

          assignments1         = load (assignments1Filename);

          assignments2         = load (assignments2Filename);

          assignments2         = compareAssignments1And2(assignments1,assignments2,modelIndex,modelIndex+1);
          if (~isempty(assignments2))
              compareAssignments1And2(assignments2,assignments1,modelIndex+1,modelIndex);
          end
          fprintf(1, '-----------------------------------------------------------\n');
          fprintf(1, '-----------------------------------------------------------\n');
      end
    end
end


%%
%is a directional comparison
function assignments2Pruned = compareAssignments1And2(assignments1,assignments2,model1Index,model2Index,fid)

for index = 1:size(assignments1)
	
  peakIndex    = assignments1(index,1);
  residueIndex = assignments1(index,2);
          
  indexInSecondAssignment = find(assignments2(:,1) == peakIndex);
  
  if (isempty(indexInSecondAssignment))
    if (model2Index < model1Index)
        fprintf (1, '%d %d %d -1 %d %d \n',model2Index,model1Index,peakIndex,peakIndex,residueIndex);
    else
        fprintf (1, '%d %d %d %d %d -1\n',model1Index,model2Index,peakIndex,residueIndex,peakIndex);
    end
  else
      
    residue2Index = assignments2(indexInSecondAssignment,2);
    if (residue2Index ~= residueIndex)
      if (model2Index < model1Index)
	fprintf (1, '%d %d %d %d %d %d\n',model2Index,model1Index,peakIndex,residue2Index,peakIndex,residueIndex);
      else
	fprintf (1, '%d %d %d %d %d %d\n',model1Index,model2Index,peakIndex,residueIndex,peakIndex,residue2Index);
      end
    end
    assignments2(indexInSecondAssignment,:) = [-1 -1];
  end
 
end

assignments2Pruned = removeComparedElements(assignments2);


%%
function assignments2Pruned = removeComparedElements(assignments2)
assignments2Pruned = zeros(1,2);
counter            = 1;

for i = 1:size(assignments2)
  if ((assignments2(i,1) == -1) && (assignments2(i,2) == -1))
  else
      assignments2Pruned(counter,1) = assignments2(i,1);
      assignments2Pruned(counter,2) = assignments2(i,2);
      counter = counter + 1;
  end
end     

if (counter == 1)
    assignments2Pruned = [];
end
