function analyzeNormalModeHD_Results

dbstop if error;

maxAssignmentAccuracy       = 0;
minAssignmentAccuracy       = 100;
maxHD_Score                 = -1000;
bestHD_Score_Accuracy       = 0;


numTemplates   = 4;
dirNames       = cell(numTemplates,1);

%dirNames{1}    = '1HEZ';
%dirNames{1}    = '1JML';
dirNames{1}    = '1EF1';
dirNames{2}    = '1VCB';
dirNames{3}    = '1H8C';
%dirNames{1}     = '1EF1';
%dirNames{1}   = '1DK8';
dirNames{4}   = '1RFA';

origAccuracies = [];

for dirIndex = 1:numTemplates
  
  cd (dirNames{dirIndex});
  cd ('Bidirectional');

  inputFilename         = 'Results/Mode.7.8.HD_vs_NMA.txt.startModel1.endModel121';
  
  H                     = load (inputFilename);
  maxAssignmentAccuracy = max(maxAssignmentAccuracy,max(H(:,3)));
  minAssignmentAccuracy = min(min(H(:,3)),minAssignmentAccuracy)

  
  [maxHD_Score_uniqueMode index_uniqueMode] = max(H(:,2));
  
  if (maxHD_Score_uniqueMode > maxHD_Score)
    maxHD_Score           = maxHD_Score_uniqueMode
    bestHD_Score_Accuracy = H(index_uniqueMode,3)
  end
  
  cd ('../');
  
  for modeIndex   = 7:11
    
    inputFilename = sprintf('Results/Mode%d_HD_vs_NMA.txt', modeIndex);
    
    H             = load (inputFilename);

    
    maxAssignmentAccuracy   = max(maxAssignmentAccuracy, max(H(:,3)))
    minAssignmentAccuracy   = min(minAssignmentAccuracy, min(H(:,3)))
    
    [maxHD_Score_uniqueMode index_uniqueMode] = max(H(:,2));
    
    if (maxHD_Score_uniqueMode > maxHD_Score)
      maxHD_Score           = maxHD_Score_uniqueMode
      bestHD_Score_Accuracy = H(index_uniqueMode,3)
    end
    
  end
  
  origAccuracies            = [origAccuracies H(6,3)];
  
  cd ('../');
end

fprintf(1, 'minAccuracy = %f maxAccuracy = %f bestHDScoreAccuracy=%f\n',minAssignmentAccuracy,maxAssignmentAccuracy,bestHD_Score_Accuracy);
origAccuracies
