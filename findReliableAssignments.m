function [resonanceAssignmentAccuracy, confidentResonanceAssignmentAccuracy, numConfidentPeaks, ...
	  numCorrectConfidentPeaks] = findReliableAssignments(includeUnidirectionalMode,includeBidirectionalMode)

if (nargin == 0)
  includeUnidirectionalMode = 1;
  includeBidirectionalMode  = 1;
elseif (nargin == 1)
  assert (includeUnidirectionalMode == 1);
  includeBidirectionalMode  = 0;
end
load ('~/Workdir/allPeaksAfterInitialAssignments-WithMBShiftAndRDC_Score-PartiallyCorrectRDC_Tensor-higherRDC_ScoreThresholds.mat');

printAssignmentAccuracy = 0;
csCoefficient = 0.19; rdcCoefficient = 0.03; noeCoefficient = 0.78;
numModelsPerMode             = 11;
highestNormalModeIndex       = 11;
totalNumModels               = 121+55;
numBidirectionalModels       = 121;

confidenceThreshold          = 0.5;
%numTemplates                 = 4; 
numTemplates                 = 1; 
%dirNames                     = cell(numTemplates,1);
assignments                  = cell(totalNumModels,1);
numModels                    = 0;

relNumResidues               = size(MASTER,2);

%dirNames{1} = '1HEZ';
%dirNames{1} = '1H8C';
%dirNames{2} = '1RFA';
%dirNames{3} = '1VCB';
%dirNames{4} = '1EF1';


%cd (dirNames{1});

%reading ideal assignments ...

load idealAssignments.txt
[numPeaks,unusedVar]      = size(idealAssignments);

numResidues               = max(idealAssignments(:,2));
%keyboard
idealAssignmentMatrix     = zeros(numPeaks, numResidues);



for peakIndex = 1:numPeaks
  idealAssignmentMatrix(peakIndex,idealAssignments(peakIndex,2)) = 1;
end

%cd ..;
assignmentAccuracies = [];MB_Scores  = []; BayesianScores       = ...
    []; marsScores= [];
%dbstop if error;

if (includeBidirectionalMode)

  for dirIndex = 1:numTemplates

%    cd (dirNames{dirIndex});
    
    for modelIndex = 1:numBidirectionalModels

assignmentsFilename          = sprintf('Bidirectional/Assignments/assignmentsForMode7.8.Model%d.txt',modelIndex);
%      assignmentsFilename       = sprintf('assignmentsForMode78_Model%d.txt', modelIndex);

      thisMASTER                = zeros(numPeaks,relNumResidues);
      
      if (~checkIfFileExists(assignmentsFilename))
	fprintf(1, 'couldnt read %s\n',assignmentsFilename);
	keyboard
	continue;
      end
      
      numModels = numModels + 1;
      
      assignments{numModels}         = load (assignmentsFilename);
      
      for peakIndex = 1:numPeaks
	
	correspondingResidueIndex    = assignments{numModels}(peakIndex,2);
	
	relResidueIndex              = find(RESNUMS == ...
					    correspondingResidueIndex);
	
	assert (~isempty(relResidueIndex));
	
	if (correspondingResidueIndex == -1)
	  %unassigned peak
	  continue;
	end
	
	thisMASTER(peakIndex,relResidueIndex) ...
	    = 1;
      end
      
      assignmentAccuracy   = computeAssignmentAccuracy(peakIDs, ...
				RESNUMS,thisMASTER, printAssignmentAccuracy);
      
      MB_Score             = computeMB_Score         (thisMASTER, MB_ShiftScore,MB_RDC_Score,numCS,numRDC,HSQCDATA,ALLDISTS,NTH,csCoefficient,rdcCoefficient,noeCoefficient);
      
      BayesianScore        = computeBayesianScore    (thisMASTER);
      
      marsScore            = computeMarsScore(thisMASTER);
      
      marsScores           = [marsScores marsScore];
      
      BayesianScores       = [BayesianScores BayesianScore];
      
      assignmentAccuracies = [assignmentAccuracies assignmentAccuracy];
      
      MB_Scores            = [MB_Scores MB_Score];
      
    end
%    cd ..;
  end
end

if (includeUnidirectionalMode)
  for dirIndex = 1:numTemplates
 %   cd (dirNames{dirIndex});
    for modeIndex = 7:highestNormalModeIndex
      for modelIndex = 1:numModelsPerMode
	thisMASTER                = zeros(numPeaks,relNumResidues);
	assignmentsFilename          =  sprintf('Assignments/assignmentsForMode%d_Model%d.txt',	modeIndex, modelIndex);
	%assignmentsFilename          = sprintf('assignmentsForMode%d_Model%d.txt', modeIndex, modelIndex);
	
	if (~checkIfFileExists(assignmentsFilename))
	  fprintf(1, 'couldnt read %s\n',assignmentsFilename);
	  keyboard
	  continue;
	end

	numModels = numModels + 1;
	
	assignments{numModels}         = load (assignmentsFilename);
	
	for peakIndex = 1:numPeaks
	  
	  correspondingResidueIndex  = assignments{numModels}(peakIndex,2);
	  
	  relResidueIndex              = find(RESNUMS == ...
					      correspondingResidueIndex);
	  
	  assert (~isempty(relResidueIndex));
	  
	  
	  if (correspondingResidueIndex == -1)
	    %unassigned peak
	    continue;
	  end
	  
	  thisMASTER(peakIndex,relResidueIndex) ...
	      =  1;
	end
	
	assignmentAccuracy   = computeAssignmentAccuracy(peakIDs, ...
							 RESNUMS,thisMASTER, printAssignmentAccuracy);
	
	MB_Score             = computeMB_Score         (thisMASTER, MB_ShiftScore,MB_RDC_Score,numCS,numRDC,HSQCDATA,ALLDISTS,NTH,csCoefficient,rdcCoefficient,noeCoefficient);
	
	BayesianScore        = computeBayesianScore    (thisMASTER);
	
	assignmentAccuracies = [assignmentAccuracies assignmentAccuracy];
	
	MB_Scores            = [MB_Scores MB_Score];
	
	BayesianScores       = [BayesianScores BayesianScore];

	marsScore            = computeMarsScore(thisMASTER);
	
	marsScores           = [marsScores marsScore];
	
	
      end
    end
%    cd ..;
  end
end

corrCoefMatrix = corrcoef(marsScores, assignmentAccuracies);
figure; 
plot(marsScores, assignmentAccuracies,'*');
xlabel('MarsScores')
ylabel('Assignment Accuracies');
fprintf(1, 'the corr coef is %f\n',corrCoefMatrix(1,2));

keyboard


corrCoefMatrix = corrcoef(BayesianScores, assignmentAccuracies);
figure; 
plot(BayesianScores, assignmentAccuracies,'*');
xlabel('BayesianScores')
ylabel('Assignment Accuracies');
fprintf(1, 'the corr coef is %f\n',corrCoefMatrix(1,2));

keyboard


corrCoefMatrix = corrcoef(MB_Scores, assignmentAccuracies);
figure; 
plot(MB_Scores, assignmentAccuracies,'*');
xlabel('MB_Scores')
ylabel('Assignment Accuracies');
fprintf(1, 'the corr coef is %f\n',corrCoefMatrix(1,2));

keyboard

numModels

correctness = zeros(numPeaks,1);
concurRatio = correctness;

figure;

correctPeaks         = [];
correctConfidences   = [];
incorrectPeaks       = [];
incorrectConfidences = [];

allConfidences       = [];

assignments          = [];

assert (numResidues >= numPeaks);

squareAssignmentMatrix                  = zeros(numResidues,numResidues);
squareAssignmentMatrix(1:numPeaks,:)    = assignmentMatrix;


 
 
 [MBM_Assignments,T]      = hungarian(-squareAssignmentMatrix');
 
 check = 0;
 for i = 1:numPeaks
   check = check + squareAssignmentMatrix(i,MBM_Assignments(i));
 end
 
 assert (T == -check);
 for i = numPeaks+1:numResidues
   check = check + squareAssignmentMatrix(i,MBM_Assignments(i));
 end

 assert (T == -check);

len_MBM_Assignments = length(MBM_Assignments);
newFormat_MBM_Assignments = zeros(len_MBM_Assignments,2);
for i = 1:len_MBM_Assignments
  newFormat_MBM_Assignments(i,1) = i; newFormat_MBM_Assignments(i,2) ...
      = MBM_Assignments(i);
end
 
 
 
for peakIndex   = 1:numPeaks

   residueIndex               = MBM_Assignments(peakIndex);

%   residueIndex               = SM_Assignments(peakIndex,2);
  
%   keyboard
   
   numConcurringAssignments   = assignmentMatrix(peakIndex,residueIndex);
   
   assignments                = [assignments residueIndex];

   peaksAssignedToThisResidue = find(assignments == ...
				     residueIndex);
  
   numPeaksAssignedToThisResidue = length(peaksAssignedToThisResidue);
   
   if (numPeaksAssignedToThisResidue > 1)
     fprintf(1, 'WARNING: the following peaks were assigned to this_residue:\n');
     
     for incorrectlyAssignedPeakIndex = 1:numPeaksAssignedToThisResidue
       fprintf(1, '%d ',peaksAssignedToThisResidue(incorrectlyAssignedPeakIndex));
     end
     
     fprintf(1, '\n');
     
     keyboard
   end

  

   
   
   correctResidueIndex                     = find(idealAssignmentMatrix(peakIndex,:));
   %note: this is an absolute index.
   
   assert (~isempty(correctResidueIndex));
  
   concurRatio(peakIndex) = numConcurringAssignments/numModels;
   allConfidences     = [allConfidences concurRatio(peakIndex)];
   
   if (residueIndex == correctResidueIndex)
     fprintf(1, 'numConcurringAssignments = %d, peak %d CORRECT assignment\n',numConcurringAssignments,peakIndex);
     fprintf(1, 'correct residue = %d\n', correctResidueIndex);
     correctness(peakIndex) = 1;
     correctPeaks       = [correctPeaks peakIndex];
     correctConfidences = [correctConfidences concurRatio(peakIndex)];
     hold on
   else
     
     fprintf(1, 'numConcurringAssignments = %d, peak %d WRONG assignment\n' ,numConcurringAssignments,peakIndex);
     fprintf(1, 'correct residue = %d, numConcurringAssignments_for_that=%d\n', correctResidueIndex,assignmentMatrix(peakIndex,correctResidueIndex));
     
     correctness(peakIndex) = 0;
     incorrectPeaks = [incorrectPeaks peakIndex];
     incorrectConfidences = [incorrectConfidences concurRatio(peakIndex)];
     hold on
   end      
end



[sortedConfidences, sortedPeakIndices] = sort(allConfidences);

newCorrectPeaks         = [];
newIncorrectPeaks       = [];
newCorrectConfidences   = [];
newIncorrectConfidences = [];
 
for newPeakIndex = 1:numPeaks
  oldPeakIndex = sortedPeakIndices(newPeakIndex);
  isCorrect    = correctness(oldPeakIndex);
  
  if (isCorrect)
    newCorrectPeaks         = [newCorrectPeaks newPeakIndex];
    newCorrectConfidences   = [newCorrectConfidences sortedConfidences(newPeakIndex)];
  else
    newIncorrectPeaks       = [newIncorrectPeaks newPeakIndex];
    newIncorrectConfidences = [newIncorrectConfidences sortedConfidences(newPeakIndex)];
  end
  
end



fid = fopen('confidencesAndCorrectness.txt','w');
fprintf(1, 'check out confidencesAndCorrectness.txt\n');

for peakIndex = 1:length(allConfidences)
  if (correctness(peakIndex))
    fprintf(fid, '%d %d %d\n', peakIndex, allConfidences(peakIndex), 1);
  else
    fprintf(fid, '%d %d %d\n', peakIndex, allConfidences(peakIndex), 0);
  end
end

fclose(fid);



resonanceAssignmentAccuracy = sum(correctness)/numPeaks


confidentAssignmentsPeakIndices      = find(concurRatio > ...
					    confidenceThreshold);

numConfidentPeaks              = length(confidentAssignmentsPeakIndices)

numCorrectConfidentPeaks       = sum(correctness(confidentAssignmentsPeakIndices))
numIncorrectConfidentAssignments     = numConfidentPeaks - numCorrectConfidentPeaks;

confidentResonanceAssignmentAccuracy = numCorrectConfidentPeaks/numConfidentPeaks






confidentAssignmentsResidueIndices   = ...
    zeros(numConfidentPeaks,1);

for i = 1:numConfidentPeaks
  peakIndex                             = confidentAssignmentsPeakIndices(i);
  [numConcurringModels,residueIndex]    = max(assignmentMatrix(peakIndex,:));
  confidentAssignmentsResidueIndices(i) = residueIndex;
end



confidentAssignmentsFilename         = 'confidentAssignments.txt';
fid                                  = fopen(confidentAssignmentsFilename,'w');
fprintf    (1,'check out %s\n',confidentAssignmentsFilename);
for i = 1:numConfidentPeaks
  fprintf    (fid,'%d %d\n',confidentAssignmentsPeakIndices(i), ...
	      confidentAssignmentsResidueIndices(i));
end
fclose     (fid);


plot(newCorrectPeaks, newCorrectConfidences, 'o', newIncorrectPeaks, ...
     newIncorrectConfidences,'rx','LineWidth',2);


legend ('correct assignments', 'incorrect','Location','Best');
xlabel('Peak index')
ylabel('confidence');

plot(0:numPeaks,confidenceThreshold)

titleString = sprintf('ratio of the models that agree on the assignment, blue:correct, red: incorrect,r.a.a=%.2f,r.a.a.for.confidentAssignments=%.2f\nnumConfidentPeaks=%d,numCorrectConfidentPeaks=%d',resonanceAssignmentAccuracy,confidentResonanceAssignmentAccuracy,numConfidentPeaks,numCorrectConfidentPeaks);

title('1RFA, confidences');

print -djpeg confidences_sorted.jpg
print -depsc confidences_sorted.eps


figure;
plot(correctPeaks, correctConfidences, 'o', incorrectPeaks, ...
     incorrectConfidences,'rx','LineWidth',2);
legend ('correct assignments', 'incorrect','Location','Best');
xlabel('Peak index')
ylabel('Percentage of models that agree on this assignment');
hold on

plot(0:numPeaks,confidenceThreshold)

titleString = sprintf('ratio of the models that agree on the assignment, blue:correct, red: incorrect,r.a.a=%.2f,r.a.a.for.confidentAssignments=%.2f\nnumConfidentPeaks=%d,numCorrectConfidentPeaks=%d',resonanceAssignmentAccuracy,confidentResonanceAssignmentAccuracy,numConfidentPeaks,numCorrectConfidentPeaks);

title(titleString);

legend('correct assignment','incorrect','Location','Best');


ylabel('confidence');
print -djpeg consensusAssignments.jpg
print -depsc consensusAssignments.eps







exit;