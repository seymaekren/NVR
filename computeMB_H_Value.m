function  [newH_Value,CS_RDC_Component] = computeMB_H_Value(newAssigntable, MASTER, HSQCDATA,...
					 unassignedPeakIndices, unassignedResidueIndices,...
					 ALLDISTS, NTH, ...
					 MB_ShiftScore, ...
					 MB_RdcScore, numCS, numRDC, ...
					 csCoefficient, rdcCoefficient, ...
					 noeCoefficient, numHN_NOES);

newH_Value = 0; 

for relPeakIndex = 1:size(newAssigntable,1)

  relResidueIndices = find(newAssigntable(relPeakIndex,:));

  peakIndex         = unassignedPeakIndices(relPeakIndex);
  
  %the peak has at least one residue it can be assigned to.
  
  assert (~isempty(relResidueIndices));

  minScore = 0;
   
  for j = 1:length(relResidueIndices)

    relResidueIndex  = relResidueIndices(j);
  
    residueIndex     =  unassignedResidueIndices(relResidueIndex);

    score            = - csCoefficient*MB_ShiftScore(peakIndex,residueIndex) ...
	- rdcCoefficient * MB_RdcScore(peakIndex,residueIndex);

    if (j == 1) | (score < minScore)
      minScore = score;
    end

%    fprintf(1, 'score of assigning peak #%d to residue#%d is %f\n',peakIndex,residueIndex,score);
    
  end
  
  newH_Value = newH_Value + minScore;
  
end

CS_RDC_Component = newH_Value;

MASTER(unassignedPeakIndices,unassignedResidueIndices) = newAssigntable;

noeScore = 0;
if (noeCoefficient > 0)
  noeScore    = computeMB_NoeScore(HSQCDATA, ALLDISTS, NTH,MASTER);
end

newH_Value  = newH_Value - noeCoefficient*noeScore;
%newH_Value  = newH_Value / normalizingConstant;

% $$$ [indexOfNoesThatMatchHN, indexOfNoesThatMatchH2]=     computeNoeHSQC_Correspondence(HSQCDATA); 
% $$$ %make noe1 and noe2 static
% $$$ assert (size(indexOfNoesThatMatchHN,1) == size(MASTER,1));
% $$$ assert (size(indexOfNoesThatMatchH2,1) == size(MASTER,1));
% $$$ 
% $$$ %keyboard
% $$$ 
% $$$ for peak1Index = 1:size(MASTER,1)
% $$$ 
% $$$   residue1Index        = find(MASTER(peak1Index,:));
% $$$ 
% $$$   if (isempty(residue1Index)) 
% $$$     relPeak1Index      = find(unassignedPeakIndices == peak1Index);
% $$$     if (isempty(relPeak1Index))%this subproblem (i.e. the assignment
% $$$                                %of this peak) is not being solved
% $$$                                %right now.
% $$$       continue;
% $$$     end
% $$$   end
% $$$   
% $$$   
% $$$   for peak2Index = 1:size(MASTER,1)
% $$$   
% $$$     residue2Index        = find(MASTER(peak2Index,:));
% $$$     
% $$$     if (~isempty(residue1Index)) & (~isempty(residue2Index))
% $$$       %this cost must already be counted in the score.
% $$$       continue;
% $$$     end
% $$$ 
% $$$     if (isempty(residue2Index)) 
% $$$       relPeak2Index      = find(unassignedPeakIndices == peak2Index);
% $$$       if (isempty(relPeak2Index))%this subproblem (i.e. the assignment
% $$$ 				 %of this peak) is not being solved
% $$$ 				 %right now
% $$$ 	continue;
% $$$       end
% $$$     end
% $$$     
% $$$     numNOEsBetweenPeaks  = findNOEsBetweenThesePeaks(indexOfNoesThatMatchHN,indexOfNoesThatMatchH2,peak1Index,peak2Index);
% $$$ 
% $$$ %    keyboard
% $$$     
% $$$     if (numNOEsBetweenPeaks == 0)
% $$$       continue;
% $$$     end
% $$$     
% $$$     
% $$$     if (isempty(residue1Index))
% $$$       assert (~isempty(relPeak1Index));
% $$$       relResidue1Indices = find(newAssigntable(relPeak1Index,:));
% $$$       residue1Indices    = unassignedResidueIndices(relResidue1Indices);
% $$$     else
% $$$       residue1Indices    = residue1Index;
% $$$     end
% $$$     
% $$$ 
% $$$     
% $$$     if (isempty(residue2Index))
% $$$ %      relPeak2Index      = find(unassignedPeakIndices ==
% $$$ %      peak2Index); This has been done above.
% $$$       assert (~isempty(relPeak2Index));
% $$$       relResidue2Indices = find(newAssigntable(relPeak2Index,:));
% $$$       residue2Indices    = unassignedResidueIndices(relResidue2Indices);
% $$$     else
% $$$       residue2Indices    = residue2Index;
% $$$     end
% $$$     
% $$$     distances            = ALLDISTS(residue1Indices, residue2Indices);
% $$$     minDistance          = min(min(distances));
% $$$     if (minDistance < NTH)
% $$$       fprintf(1, 'peak1 = %d peak2 = %d numNOEs between them = %d contribute.\n',peak1Index,peak2Index,numNOEsBetweenPeaks);
% $$$       newH_Value         = newH_Value - noeCoefficient * numNOEsBetweenPeaks;
% $$$     end
% $$$   end
% $$$ end

newH_Value = newH_Value / (csCoefficient*numCS +  rdcCoefficient*numRDC ...
			   + noeCoefficient*numHN_NOES);
