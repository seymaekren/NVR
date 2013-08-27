load ('allPeaksAfterInitialAssignments-WithMBShiftAndRDC_Score-PartiallyCorrectRDC_Tensor-higherRDC_ScoreThresholds.mat');

%keyboard

numAdditions = 0; csWeight = 0; rdcWeight = 0; noeWeight = 0;

[numPeaks,numResidues] = size(MASTER);

myAlmostCompleteMASTER = MASTER*0;

for peakIndex = 1:numPeaks
  myAlmostCompleteMASTER(peakIndex,peakIndex) =  1;
end
[correctAssignmentNoeScore,numHN_NOEs] = computeMB_NoeScore(HSQCDATA,ALLDISTS,NTH,myAlmostCompleteMASTER);
correctAssignmentNoeScore              = correctAssignmentNoeScore;%/numHN_NOEs;


%normalizedMB_Noe_Score = zeros(numPeaks,numResidues);
if (isempty(MB_Noe_Score))
  MB_Noe_Score = zeros(numPeaks,numResidues);
  for peakIndex = 1:numPeaks
    
    fprintf(1, '.');
    if (peakIndex > 1)
      myAlmostCompleteMASTER(peakIndex-1,:)           = 0;
      myAlmostCompleteMASTER(peakIndex-1,peakIndex-1) = 1;
    end
    
    for residueIndex = 1:numResidues
      if (residueIndex == peakIndex)
	MB_NoeScore(peakIndex,peakIndex)             = correctAssignmentNoeScore;
	continue;
      end
      
      
      myAlmostCompleteMASTER(peakIndex,:)            = 0;
      myAlmostCompleteMASTER(peakIndex,residueIndex) = 1;
      
      newNoeScore                                    = computeMB_NoeScore(HSQCDATA,ALLDISTS,NTH,myAlmostCompleteMASTER);
      
      %    normalizedMB_NoeScore(peakIndex,residueIndex) = newNoeScore/numHN_NOEs;
      MB_Noe_Score(peakIndex,residueIndex)           = newNoeScore;
      
    end
  end
  fprintf(1, '\n');
  
  save pairwiseMB_NoeScore.mat
else
  load pairwiseMB_NoeScore.mat
end
%csCoefficientList   = [0:0.1:1];
%rdcCoefficientList  = [0:0.1:1];  
%noeCoefficientList  = [0:0.1:1];

%csCoefficientList   = [0:0.1:1];
%rdcCoefficientList  = [0:0.1:1];  
%noeCoefficientList  = [0:0.1:1];

csCoefficientList   = [0:0.25:1];
rdcCoefficientList  = [0:0.25:1];  
noeCoefficientList  = [0:0.25:1];

coefficientIsBest  = zeros(length(csCoefficientList),length(rdcCoefficientList),length(noeCoefficientList));

%normalizedMB_Shift_Score = MB_Shift_Score/numCS;
%normalizedMB_RDC_Score   = MB_RDC_Score/numRDC;

numCases = 0;

for csCoefficientIndex = 1:length(csCoefficientList)
  fprintf(1, '.');
  csCoefficient = csCoefficientList(csCoefficientIndex);
  for rdcCoefficientIndex = 1:length(rdcCoefficientList)
    rdcCoefficient = rdcCoefficientList(rdcCoefficientIndex);
    for noeCoefficientIndex = 1:length(noeCoefficientList)
      noeCoefficient = noeCoefficientList(noeCoefficientIndex);
      for peakIndex = 1:numPeaks
	for residueIndex = 1:numResidues
	  MB_Score(peakIndex,residueIndex) = csCoefficient* ...
	      MB_ShiftScore(peakIndex,residueIndex) + ...
	      rdcCoefficient*MB_RDC_Score(peakIndex,residueIndex) ...
	      + noeCoefficient * MB_Noe_Score(peakIndex, ...
					     residueIndex);
	  normalizingConstant = csCoefficient*numCS+rdcCoefficient*numRDC+noeCoefficient*numHN_NOEs;
	  MB_Score(peakIndex,residueIndex) = MB_Score(peakIndex,residueIndex)/normalizingConstant;
	end
	
	if (MB_Score(peakIndex,peakIndex) == ...
	    max(MB_Score(peakIndex,:)))
	  
%	  fprintf(1, 'peak %d correctly assigned with coeffs\n',peakIndex);
%	  fprintf(1, '%f %f %f\n',csCoefficientList(csCoefficientIndex),rdcCoefficientList(rdcCoefficientIndex),noeCoefficientList(noeCoefficientIndex));
	  
%	  keyboard
	  coefficientIsBest(csCoefficientIndex, ...
			    rdcCoefficientIndex,noeCoefficientIndex) ...
	      = coefficientIsBest(csCoefficientIndex, ...
				  rdcCoefficientIndex,noeCoefficientIndex) + 1;
	end
	numCases = numCases + 1;
      end
    end
  end
end

maxValue                                                  = max(max(max(coefficientIsBest)));
linearIndexOfBestCoefficient                              = find(coefficientIsBest == maxValue);
[bestCsCoefficient,bestRdcCoefficient,bestNoeCoefficient] = ind2sub(size(coefficientIsBest),linearIndexOfBestCoefficient);

fprintf(1, 'percentage of time a coeff. combination ranked correctly');
fprintf(1, ' is %f\n',sum(sum(sum(coefficientIsBest)))/numCases);

%[bestCsCoefficient,bestRdcCoefficient,bestNoeCoefficient] = ...
%    find(coefficientIsBest == maxValue);
for coefficientIndex = 1:length(bestCsCoefficient)
  fprintf(1, 'best cs coefficient %f best rdc Coefficient %f ',csCoefficientList(bestCsCoefficient(coefficientIndex)),rdcCoefficientList(bestRdcCoefficient(coefficientIndex)));
  fprintf(1, 'best noe coefficient %f\n',noeCoefficientList(bestNoeCoefficient(coefficientIndex)));
  fprintf(1, 'The ratio of cases in which this combination' );
  fprintf(1, ' ranked the correct assignment first=%f\n',maxValue/numCases);
end



figure; 

