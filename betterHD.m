dbstop if error
dbstop if warning


if (refineWithRDCs)
  fprintf(1, 'refining with RDCs...\n');
else
  fprintf(1, 'not refining with RDCs...\n');
end

if (~exist('ROWIN'))

  if (useCH_RDCs)
    [ROWIN, COLIN, ASSIGNTABLE, MASTER, HDE, TP, CP, SXCP, SSCP, RDC1, RDC2, NTH, ...
     incorrectRP1, incorrectRP2, incorrectS1, incorrectS2, marsShiftScore, ...
     MB_ShiftScore, numCS] = initialize(peaks,rdcs,...
					HDEXCHANGE, ...
					peakIDs, NOES, ...
					VECTORS_NH,TYPES, ...
					RESNUMS,SSTRUCT, HBOND, ...
					ALLDISTS,IALLDISTS, ...
					SHIFTS_Filename, SHIFTX_Filename, ...
					useCH_RDCs,useHD_Routines, useTOCSY, ...
					truncateProbabilities);

  elseif ((refineWithRDCs == 1) & (b_runningMBP == 1))
    [ROWIN, COLIN, ASSIGNTABLE, MASTER, HDE, TP, CP, SXCP, SSCP, RDC1, RDC2, NTH, ...
     incorrectRP1, incorrectRP2, incorrectS1, incorrectS2, marsShiftScore, ...
     MB_ShiftScore, numCS] = initialize(peaks,rdcs,...
					HDEXCHANGE, ...
					peakIDs, NOES, ...
					VECTORS_NH,TYPES, ...
					RESNUMS,SSTRUCT, HBOND, ...
					ALLDISTS,IALLDISTS, ...
					SHIFTS_Filename, SHIFTX_Filename, ...
					useCH_RDCs,useHD_Routines, useTOCSY, ...
					truncateProbabilities);
  else
  
    [ROWIN, COLIN, ASSIGNTABLE, MASTER, HDE, TP, CP, SXCP, SSCP, RDC1, RDC2, NTH, RP1, RP2, S1, S2, marsShiftScore, MB_ShiftScore, numCS] = initialize(peaks,rdcs,...
						  HDEXCHANGE, ...
						  peakIDs, NOES, ...
						  VECTORS,TYPES, ...
						  RESNUMS,SSTRUCT, HBOND, ...
						  ALLDISTS,IALLDISTS, ...
						  SHIFTS_Filename, ...
						  SHIFTX_Filename, ...
						  useCH_RDCs,useHD_Routines, useTOCSY, truncateProbabilities);
    
  end
 
  b_reduceSizeOfTheProblem = 0;

  if (b_reduceSizeOfTheProblem)  

    reduceSizeOfTheProblem();
    
  end
  
  if (b_printOptimizationFiles) | (b_printIndividualScores)
    
    myMASTER                  = MASTER*0;
    numPeaks                  = size(myMASTER,1);

    if (refineWithRDCs)
      myMASTER_WithExtraColumns =     load (assignmentMatrixFilename);
      fprintf(1, 'reading %s...\n', assignmentMatrixFilename);
      fprintf(1, 'the matrix read has %d rows and %d columns.\n', ...
	      size(myMASTER_WithExtraColumns,1), size(myMASTER_WithExtraColumns,2));
      myMASTER                  =     myMASTER_WithExtraColumns(:,1:size(myMASTER_WithExtraColumns,2)-3);
      fprintf(1, 'the derived matrix has %d columns.\n', size(myMASTER,2));
%      fprintf(1, 'KEYBOARD COMMENTED OUT FOR THE AUTOMATED TEST.\n');
      keyboard
      
      numCorrectAssignmentsInTheFileRead = 0;

%      fprintf(1, 'WARNING, SETTING THE READ MATRIX TO BE 1 ON DIAGONAL.\n');
	
      for peakIndex = 1:numPeaks
%	myMASTER(peakIndex,:) = 0;
%	myMASTER(peakIndex,peakIndex) = 1;
	if (myMASTER(peakIndex,peakIndex) == 1)
	  numCorrectAssignmentsInTheFileRead = ...
	      numCorrectAssignmentsInTheFileRead + 1;
	end
      end
      fprintf(1, 'assignment accuracy in the inputted file = %f\n',numCorrectAssignmentsInTheFileRead/size(myMASTER,1));

      if (useCH_RDCs == 1)
	S1                               = updateTen_CH(myMASTER, ...
							NH_RDCS,CH_RDCS, VECTORS_NH, VECTORS_CH);
	RP1                              = NVR_RDC2PROB(ASSIGNTABLE,NH_RDCS,VECTORS_NH,S1,ROWIN,COLIN);
	RP2                              = NVR_RDC2PROB(ASSIGNTABLE,CH_RDCS,VECTORS_CH,S1,ROWIN,COLIN);

	if (useFourRDCsPerResidue)
	  S2                               = updateTen(myMASTER,RDC1,VECTORS_NH);
	  S3                               = updateTen(myMASTER,RDC2,VECTORS_NH);
	  RP3                              = NVR_RDC2PROB(ASSIGNTABLE,RDC1,VECTORS_NH,S2,ROWIN,COLIN);
	  RP4                              = NVR_RDC2PROB(ASSIGNTABLE,RDC2,VECTORS_NH,S3,ROWIN,COLIN);
	end
      end

      if (useCH_RDCs == 0)

	if (b_runningMBP == 1)
	  assert (refineWithRDCs == 1);
%	  S1                               = updateTen_CCa_NH_NC(myMASTER, ...
%						  CCa_RDCS, ...
%						  NH_RDCS, NC_RDCS, ...
%						  VECTORS_CCa, VECTORS_NH, VECTORS_NC);
          S1                               = updateTen(myMASTER,CCa_RDCS,VECTORS_CCa)
%	  S2                               = updateTen(myMASTER,NH_RDCS,VECTORS_NH)
	  S3                               = updateTen(myMASTER,NC_RDCS,VECTORS_NC)
          keyboard
	  RP1                              = NVR_RDC2PROB(ASSIGNTABLE,CCa_RDCS,VECTORS_CCa,S1,ROWIN,COLIN);
%          RP2                              = NVR_RDC2PROB(ASSIGNTABLE,NH_RDCS,VECTORS_NH,S2,ROWIN,COLIN);
	  RP3                              = NVR_RDC2PROB(ASSIGNTABLE,NC_RDCS,VECTORS_NC,S3,ROWIN,COLIN);
	else
	  S1                               = updateTen(myMASTER,RDC1,VECTORS);
	  S2                               = updateTen(myMASTER,RDC2,VECTORS);
	  RP1                              = NVR_RDC2PROB(ASSIGNTABLE,RDC1,VECTORS,S1,ROWIN,COLIN);
	  RP2                              = NVR_RDC2PROB(ASSIGNTABLE,RDC2,VECTORS,S2,ROWIN,COLIN);
	  %	printTen(S1, 's1Tensor.txt');
	  %	printTen(S2, 's2Tensor.txt');
	  %	keyboard
	end
      end
    else
      RP1 = []; RP2 = [];
    end



    useBayesianMatrix                = 0;
    
    if (useBayesianMatrix)

      numNvrVoters                     = 5;
      nvrVoter                         = initializeVoters(SSCP,SXCP,CP,TP,HDE, RP1, RP2);
      bayesianMatrix                   = computeBayesianMatrix(ASSIGNTABLE, ROWIN, COLIN, ...
						  nvrVoter, numNvrVoters);
      fprintf(1, 'size of bayesianMatrix is %d and %d\n', size(bayesianMatrix,1),size(bayesianMatrix,2));
    end
    
    if (b_printOptimizationFiles)
    
      if (useFourRDCsPerResidue)
      
	printOptimizationFiles_4RDC_Matrices(CP, SXCP, SSCP, TP, HDE, ...
					     RP1, RP2, RP3, RP4, NOES, ALLDISTS, COLIN, NTH);
	
      elseif ((b_runningMBP == 1) & (refineWithRDCs == 1))
%	printOptimizationFiles_3RDC_Matrices(CP, SXCP, SSCP, TP, HDE, ...
%					     RP1, RP2, RP3, NOES, ALLDISTS, COLIN, NTH);
      	printOptimizationFiles(CP, SXCP, SSCP, TP, HDE, refineWithRDCs, ...
			       RP1, RP3, NOES, ALLDISTS, COLIN, NTH, ...
			       useBayesianMatrix);
      elseif (b_printIndividualScoringMatrix)

	if (b_printCP)
	  scoringMatrixFilename = 'individualScoringMatrix_CP.txt';
	  printOptimizationFiles_IndScoringMatrix(CP, NOES, ALLDISTS, ...
						  COLIN, NTH, scoringMatrixFilename);
	elseif (b_printSXCP)
	  scoringMatrixFilename = 'individualScoringMatrix_SXCP.txt';
	  printOptimizationFiles_IndScoringMatrix(SXCP, NOES, ALLDISTS, COLIN, NTH, scoringMatrixFilename);
	elseif (b_printSSCP)
	  scoringMatrixFilename = 'individualScoringMatrix_SSCP.txt';
	  printOptimizationFiles_IndScoringMatrix(SSCP, NOES, ALLDISTS, COLIN, NTH, scoringMatrixFilename);
	elseif (b_printTP)
	  scoringMatrixFilename = 'individualScoringMatrix_TP.txt';
	  printOptimizationFiles_IndScoringMatrix(TP, NOES, ALLDISTS, COLIN, NTH, scoringMatrixFilename);
	elseif (b_printHDE)
	  scoringMatrixFilename = 'individualScoringMatrix_HDE.txt';
	  printOptimizationFiles_IndScoringMatrix(HDE, NOES, ALLDISTS, COLIN, NTH, scoringMatrixFilename);
	elseif (b_printRP1)
	  assert (refineWithRDCs == 1);
	  scoringMatrixFilename = 'individualScoringMatrix_RP1.txt';
	  printOptimizationFiles_IndScoringMatrix(RP1, NOES, ALLDISTS, COLIN, NTH, scoringMatrixFilename);
	elseif (b_printRP2)
	  assert (refineWithRDCs == 1);
	  scoringMatrixFilename = 'individualScoringMatrix_RP2.txt';
	  printOptimizationFiles_IndScoringMatrix(RP2, NOES, ALLDISTS, COLIN, NTH, scoringMatrixFilename);
	elseif (b_printCS)
	  scoringMatrixFilename = 'individualScoringMatrix_CS.txt';
	  printOptimizationFiles_IndScoringMatrix(CP, NOES, ALLDISTS, ...
						  COLIN, NTH, scoringMatrixFilename, SXCP, SSCP);
	elseif (b_printRDC)
	  assert (refineWithRDCs == 1);
	  scoringMatrixFilename = 'individualScoringMatrix_RDC.txt';
	  printOptimizationFiles_IndScoringMatrix(RP1, NOES, ALLDISTS, ...
						  COLIN, NTH, scoringMatrixFilename, RP2);
	end  
	
      else

	printOptimizationFiles(CP, SXCP, SSCP, TP, HDE, refineWithRDCs, ...
			       RP1, RP2, NOES, ALLDISTS, COLIN, NTH, ...
			       useBayesianMatrix);
      end
	
    elseif (b_printIndividualScores)
      printIndividualScores(CP, SXCP, SSCP, TP, HDE, refineWithRDCs, ...
			     RP1, RP2);
      filterScores
    end

% $$$     save ('allPeaksAssignmentEnvironment-For1UD7-forBayesianScoring-withNVR_ScoringMatrices.mat');
% $$$      
% $$$     fprintf(1, 'saved allPeaksAssignmentEnvironment-For1UD7-forBayesianScoring-withNVR_ScoringMatrices.mat\n');

  end
  
% $$$   origCP = CP; origSXCP = SXCP; origSSCP= SSCP; origTP = TP; origHDE = HDE; 
% $$$   [MASTER, ASSIGNTABLE, CP, SXCP, SSCP, RP1, RP2, TP, HDE, ROWIN, COLIN,marsRdcScore,MB_RDC_Score,numRDC]= ...
% $$$       AssignFirstFiveResidues(MASTER, ASSIGNTABLE, NOES, IALLDISTS, ...
% $$$ 			      NTH, ROWIN, COLIN, ALLDISTS, CP, SXCP, ...
% $$$ 			      SSCP, RP1, RP2, TP, HDE, S1, S2, RDC1, ...
% $$$ 			      RDC2, VECTORS);
% $$$   
% $$$   CP = origCP(ROWIN,COLIN); SXCP = origSXCP(ROWIN,COLIN); SSCP = ...
% $$$        origSSCP (ROWIN, COLIN); TP = origTP(ROWIN, COLIN); HDE = ...
% $$$        origHDE(ROWIN, COLIN); 
% $$$ %not doing something similar to the above for the RP1 and RP2 since
% $$$ %5 peaks have only been assigned now.
% $$$ 
% $$$   assert (size(RP1,1) == length(ROWIN));
% $$$   assert (size(RP1,2) == length(COLIN));
% $$$   assert (size(RP2,1) == length(ROWIN));
% $$$   assert (size(RP2,2) == length(COLIN));
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ %  save bestMarsScoreWithFullyCorrectRDCs.mat
% $$$ %  keyboard
% $$$   
% $$$   %  save allPeaksAfterInitialAssignment.mat
% $$$ %   save 20PeaksAfterInitialAssignment.mat
% $$$ %   save 11PeaksAfterInitialAssignment.mat
% $$$ % save 30PeaksAfterInitialAssignment.mat
% $$$ %  save 20PeaksAfterInitialAssignment--MoreNOEs.mat
% $$$ %  save allPeaksAfterInitialAssignments-WithMB_Scores.mat
% $$$ %save allPeaksAfterInitialAssignments-WithMB_Scores-AndStrangeSHIFTS_Matrices.mat
% $$$ %save allPeaksAfterInitialAssignments--ForUBQ.mat
% $$$ %  save allPeaksAfterInitialAssignments--For1UBI-WithFullyCorrectAlignmentTensorsForMB_RDC_Score.mat
% $$$ 
% $$$    save allPeaksAfterInitialAssignments--For1UBI-WithNVR_Matrices.mat
% $$$    fprintf(1, 'saved the environment to allPeaksAfterInitialAssignments--For1UBI-WithNVR_Matrices.mat\n');
% $$$    keyboard
end
% $$$ 
% $$$ 
% $$$ % $$$ debugNOES = NOES;
% $$$ % $$$ numPeaks  = size(NOES,1);
% $$$ % $$$ load order.m
% $$$ % $$$ 
% $$$ % $$$ for i = 1:size(ALLDISTS,1)
% $$$ % $$$   for j = 1:size(ALLDISTS,2)
% $$$ % $$$     if (ALLDISTS(i,j) < NTH)
% $$$ % $$$       residueIndex1 = order(i); residueIndex2 = order(j);
% $$$ % $$$       %fprintf(1, 'could write an NOE between residues %d and %d\n',residueIndex1,residueIndex2);
% $$$ % $$$       %assert (~isempty(peakIndex1)); assert (~isempty(peakIndex2));
% $$$ % $$$       if ((i <= numPeaks) & (j <= numPeaks))
% $$$ % $$$ 	if (NOES(i,j) == 1)
% $$$ % $$$ 	  %fprintf(1, 'set this NOE between peaks %d and %d already.\n',i,j);
% $$$ % $$$ 	  if (i < j)
% $$$ % $$$ 	    fprintf(1, 'residue %d and %d have an NOE between them.\n', ...
% $$$ % $$$ 		    residueIndex1,residueIndex2);
% $$$ % $$$ 	  end
% $$$ % $$$ 	  debugNOES (i,j) = 0;
% $$$ % $$$ 	else
% $$$ % $$$ 	  fprintf(1, 'could write an NOE between residues %d and %d\n',residueIndex1,residueIndex2);
% $$$ % $$$ 	end
% $$$ % $$$       end
% $$$ % $$$     end
% $$$ % $$$   end
% $$$ % $$$ end
% $$$ % $$$ 
% $$$ % $$$ 
% $$$ % $$$ [remainingNOEPeaks,remainingNOEResidues] = find(debugNOES);
% $$$ % $$$ assert (isempty(remainingNOEPeaks));
% $$$ % $$$ assert (isempty(remainingNOEResidues));
% $$$ % $$$ keyboard
% $$$ 
% $$$ 
% $$$ fprintf(1, 'assigned first %d residues. These are:\n', ...
% $$$ 	sum(sum(MASTER)));
% $$$ for i = 1:size(MASTER,1)
% $$$   v = find(MASTER(i,:));
% $$$   if isempty(v)
% $$$     continue;
% $$$   end
% $$$   for j = 1:length(v)
% $$$     fprintf(1, '%d %d\n', i, v(j));
% $$$   end
% $$$ end
% $$$ 
% $$$ %keyboard
% $$$ 
% $$$ 
% $$$ % $$$ if (0)
% $$$ % $$$ 
% $$$ % $$$   RP1                              = NVR_RDC2PROB(ASSIGNTABLE,RDC1,VECTORS,S1);
% $$$ % $$$   RP2                              = NVR_RDC2PROB(ASSIGNTABLE,RDC2,VECTORS,S2);
% $$$ % $$$   
% $$$ % $$$   
% $$$ % $$$   [MASTER, ASSIGNTABLE, CP, SXCP, SSCP, TP, HDE, ROWIN, COLIN]= ...
% $$$ % $$$       AssignFirstFiveResidues(MASTER, ASSIGNTABLE, NOES, IALLDISTS, ...
% $$$ % $$$ 			      NTH, ROWIN, COLIN, ALLDISTS, CP, SXCP, ...
% $$$ % $$$ 			      SSCP, RP1, RP2, TP, HDE, S1, S2, RDC1, ...
% $$$ % $$$ 			      RDC2, VECTORS);
% $$$ % $$$   
% $$$ % $$$   fprintf(1, 'assigned second set of residues. Total %d peaks assigned. These are:\n', ...
% $$$ % $$$ 	  sum(sum(MASTER)));
% $$$ % $$$   for i = 1:size(MASTER,1)
% $$$ % $$$     v = find(MASTER(i,:));
% $$$ % $$$     if isempty(v)
% $$$ % $$$       continue;
% $$$ % $$$     end
% $$$ % $$$     for j = 1:length(v)
% $$$ % $$$       fprintf(1, '%d %d\n', i, v(j));
% $$$ % $$$     end
% $$$ % $$$   end
% $$$ % $$$   
% $$$ % $$$   keyboard
% $$$ % $$$ end
% $$$ 
% $$$ 
% $$$ 
% $$$ % $$$ RP1                              = NVR_RDC2PROB(ASSIGNTABLE,RDC1,VECTORS,S1);
% $$$ % $$$ RP2                              = NVR_RDC2PROB(ASSIGNTABLE,RDC2,VECTORS,S2);
% $$$ % $$$ 
% $$$ % $$$ 
% $$$ % $$$ [MASTER, ASSIGNTABLE, CP, SXCP, SSCP, TP, HDE, ROWIN, COLIN]= ...
% $$$ % $$$       AssignFirstFiveResidues(MASTER, ASSIGNTABLE, NOES, IALLDISTS, ...
% $$$ % $$$ 			      NTH, ROWIN, COLIN, ALLDISTS, CP, SXCP, ...
% $$$ % $$$ 			      SSCP, RP1, RP2, TP, HDE, S1, S2, RDC1, ...
% $$$ % $$$ 			      RDC2, VECTORS);
% $$$ % $$$  
% $$$ % $$$ fprintf(1, 'assigned third set of residues. Total %d peaks assigned. These are:\n', ...
% $$$ % $$$ 	sum(sum(MASTER)));
% $$$ % $$$ for i = 1:size(MASTER,1)
% $$$ % $$$   v = find(MASTER(i,:));
% $$$ % $$$   if isempty(v)
% $$$ % $$$     continue;
% $$$ % $$$   end
% $$$ % $$$   for j = 1:length(v)
% $$$ % $$$     fprintf(1, '%d %d\n', i, v(j));
% $$$ % $$$   end
% $$$ % $$$ end
% $$$ % $$$ 
% $$$ % $$$ keyboard
% $$$ 
% $$$ [h_ppm n_ppm h_ppm2] = textread('dnns.txt','%f %f %f');
% $$$ 
% $$$ numHN_NOES = length(h_ppm);
% $$$ 
% $$$ csCoefficient = 0.22; rdcCoefficient = 0.05; noeCoefficient =0.73
% $$$ 
% $$$ 
% $$$ 
% $$$ [totalNumPeaks,totalNumResidues] = size(MASTER);
% $$$ numAssignedPeaks                 = 0;
% $$$ numUnassignedPeaks               = totalNumPeaks;
% $$$ 
% $$$ numNvrVoters                     = 7;
% $$$ %numVoters                        = 5;
% $$$ notCompletedAssignments          = 1;
% $$$ prevNumAssignedPeaks             = sum(sum(MASTER));
% $$$ 
% $$$ fullMASTER = MASTER*0;
% $$$ for peakIndex = 1:size(MASTER,1)
% $$$   fullMASTER(peakIndex,peakIndex) = 1;
% $$$   assert (sum(fullMASTER(peakIndex,:)) == 1);
% $$$ end
% $$$ 
% $$$ S1_full                                       = updateTen(fullMASTER, ...
% $$$ 						  RDC1, VECTORS);
% $$$ 
% $$$ S2_full                                       = updateTen(fullMASTER, ...
% $$$ 						  RDC2, VECTORS);
% $$$ 
% $$$ percentiles1 = []; percentiles2= [];
% $$$ numCorrectPerIter   = []; numTotalAssignmentsPerIter = [];
% $$$ initialMASTER = MASTER;   
% $$$ initialASSIGNTABLE = ASSIGNTABLE;
% $$$ initialROWIN = ROWIN;
% $$$ initialCOLIN = COLIN;
% $$$ initialSSCP = SSCP;
% $$$ initialSXCP = SXCP;
% $$$ initialCP = CP;
% $$$ initialTP = TP;
% $$$ initialHDE = HDE; 
% $$$ secondRun = 1;
% $$$ 
% $$$ load alignmentTensorsAfterOneRoundOfCompleteAssignments.mat
% $$$   
% $$$ %S1                               = updateTen(MASTER,RDC1,VECTORS);
% $$$ %S2                               = updateTen(MASTER,RDC2,VECTORS);
% $$$  
% $$$ %S1
% $$$ %S2
% $$$ %keyboard
% $$$ 
% $$$ 
% $$$ while (notCompletedAssignments)
% $$$ % $$$   
% $$$ % $$$   MASTER_For_BayesianMatrixComputation = MASTER*0;
% $$$ % $$$   
% $$$ % $$$   for peakIndex=1:size(MASTER,1)
% $$$ % $$$     MASTER_For_BayesianMatrixComputation(peakIndex,peakIndex) = 1;
% $$$ % $$$   end
% $$$ % $$$   
% $$$ % $$$   S1                               = updateTen(MASTER_For_BayesianMatrixComputation,RDC1,VECTORS);
% $$$ % $$$   S2                               = updateTen(MASTER_For_BayesianMatrixComputation,RDC2,VECTORS);
% $$$   
% $$$   
% $$$ 
% $$$ %here isn't S1 wrong? I mean, does it not come from a perfect
% $$$ %assignment ? It looks like it... What does this do to the
% $$$ %computation? It should compute an RP1 that is better than it
% $$$ %actually should be which is then used in BayesianScoreMatrix
% $$$ %computation. SB  
% $$$   RP1                              = NVR_RDC2PROB(ASSIGNTABLE,RDC1,VECTORS,S1,ROWIN,COLIN);
% $$$   RP2                              = NVR_RDC2PROB(ASSIGNTABLE,RDC2,VECTORS,S2,ROWIN,COLIN);
% $$$ 
% $$$   percentile1                               = ...
% $$$       NVR_COMP_TEN(S1,S1_full);
% $$$   
% $$$   percentile2                               =     NVR_COMP_TEN(S2,S2_full);
% $$$   
% $$$   fprintf(1, 'the percentile difference between the partial and full');
% $$$   fprintf(1, ' is %f and %f\n', percentile1, percentile2);
% $$$ 
% $$$   percentiles1 = [percentiles1 percentile1];
% $$$   percentiles2 = [percentiles2 percentile2];
% $$$ 
% $$$   keyboard
% $$$   
% $$$   
% $$$ 
% $$$ % $$$   fprintf(1, 'computing the Bayesian matrix using fully assigned MASTER.\n');
% $$$ % $$$   keyboard;
% $$$   
% $$$   %  marsRdcScore   = computeMarsRDC_Score(ASSIGNTABLE, RDC1, RDC2, ...
% $$$ %					VECTORS, S1, S2);
% $$$ 
% $$$ %  marsScore                        = 3.3 * marsShiftScore + marsRdcScore;
% $$$   
% $$$ %  save allPeaksAfterInitialAssignments-WithMBShiftAndRDC_Score-PartiallyCorrectRDC_Tensor-higherRDC_ScoreThresholds.mat
% $$$ 
% $$$ % $$$   fprintf(1, 'mars score computed. keyboarding.\n');
% $$$ % $$$ 
% $$$ % $$$   keyboard
% $$$ % $$$   
% $$$ % $$$   marsShiftScore = marsShiftScore(ROWIN, COLIN);
% $$$ % $$$   marsRdcScore   = marsRdcScore  (ROWIN, COLIN);
% $$$ % $$$   marsScore      = marsScore     (ROWIN, COLIN);
% $$$ % $$$   
% $$$ % $$$   LARGE_VALUE = 1E15;
% $$$ % $$$ 
% $$$ % $$$   for relPeakIndex = 1:size(marsScore,1)
% $$$ % $$$     for relResidueIndex = 1:size(marsScore,2)
% $$$ % $$$       if (marsScore(relPeakIndex,relResidueIndex) == 0)
% $$$ % $$$ 	fprintf(1, 'replacing a 0 entry.\n');
% $$$ % $$$ 	marsScore(relPeakIndex,relResidueIndex) = LARGE_VALUE;
% $$$ % $$$       end
% $$$ % $$$     end
% $$$ % $$$   end
% $$$   
% $$$   
% $$$   
% $$$   %change in variable names:
% $$$   peakIndices                      = ROWIN;
% $$$   residueIndices                   = COLIN;
% $$$ %  nvrVoter                         = cell(numVoters,1);
% $$$   nvrVoter                         = initializeVoters(SSCP,SXCP,CP,TP,HDE, RP1, RP2);
% $$$ 
% $$$ %  numVoters                        = 1;
% $$$ %  voter                            = cell(1,1);
% $$$ %  voter{1}                         = marsScore;
% $$$ %  MB_Score                         = zeros(totalNumPeaks,totalNumResidues);
% $$$ 
% $$$ %  for peakIndex = 1:totalNumPeaks
% $$$ %    for residueIndex = 1:totalNumResidues
% $$$ %      MB_Score(peakIndex,residueIndex)  = (2*MB_ShiftScore(peakIndex,residueIndex)+ MB_RDC_Score(peakIndex,residueIndex))/(2*numCS + numRDC);
% $$$ %       MB_Score(peakIndex,residueIndex)  = (MB_ShiftScore(peakIndex,residueIndex))/(numCS);
% $$$ %       voter{1}(peakIndex,residueIndex)  = exp(MB_Score(peakIndex,residueIndex));
% $$$ %       if (MB_Score(peakIndex,residueIndex) == 0)
% $$$ %	 voter{1}(peakIndex,residueIndex)  = 0;
% $$$ %       end
% $$$ %    end
% $$$  %   voter{1}(peakIndex,:) = voter{1}(peakIndex,:)/sum(voter{1}(peakIndex,:));
% $$$ %  end
% $$$ %  keyboard
% $$$ %  MB_Score = MB_Score(ROWIN,COLIN);
% $$$ %  voter{1} = voter{1}(ROWIN,COLIN);
% $$$ %  for peakIndex = 1:size(voter{1},1)
% $$$ %    voter{1}(peakIndex,:) = voter{1}(peakIndex,:)/sum(voter{1}(peakIndex,:));
% $$$ %  end
% $$$ %  voter{1} = SSCP;
% $$$ 
% $$$ %numVoters = 0;
% $$$ %voter     = cell(0,1);
% $$$ 
% $$$ 
% $$$     bayesianMatrix = computeBayesianMatrix(ASSIGNTABLE, ROWIN, COLIN, ...
% $$$   					 nvrVoter, numNvrVoters);
% $$$ % $$$  
% $$$ % $$$    
% $$$ % $$$    fullBayesianMatrix               = MASTER*0;
% $$$ % $$$    
% $$$ % $$$    fullBayesianMatrix(ROWIN, COLIN) = bayesianMatrix;
% $$$ % $$$    
% $$$ % $$$    save bayesianMatrix.mat fullBayesianMatrix
% $$$ % $$$    
% $$$ % $$$    fprintf(1, 'saved full bayesian matrix. keyboarding.\n');
% $$$ % $$$   
% $$$ % $$$   keyboard
% $$$   
% $$$   
% $$$    numVoters      = 1;
% $$$    voter          = cell(numVoters,1);
% $$$    voter{1}       = bayesianMatrix;
% $$$ 
% $$$ 
% $$$ %csCoefficient = 1; rdcCoefficient = 0; noeCoefficient =0;
% $$$ 
% $$$ %keyboard
% $$$   
% $$$   
% $$$   %  voter                            = initialize5Voters(voter,SSCP,SXCP,CP,TP,HDE);
% $$$ 
% $$$ %analyzeVoters(voter, numVoters, peakIndices, residueIndices);
% $$$ 
% $$$ %parameters of the program...
% $$$ 
% $$$ %useMBM_EM = 1
% $$$ %debugIncorrectAssignment = 0;
% $$$ 
% $$$ 
% $$$ 
% $$$ %INITIAL_MAX_NUM_RETURNED_ASSIGNMENTS = 100;
% $$$ %totalNumAssignments = 0;
% $$$ %foundMASTERs = cell(INITIAL_MAX_NUM_RETURNED_ASSIGNMENTS,1);
% $$$ %foundAnAssignment = 0;
% $$$ %minScore = 0;
% $$$ %overallMASTER = MASTER;
% $$$ %scoreSoFar = 0;
% $$$ %finalScores = zeros(INITIAL_MAX_NUM_RETURNED_ASSIGNMENTS, 1);
% $$$ %assignmentAccuracies = zeros(INITIAL_MAX_NUM_RETURNED_ASSIGNMENTS, 1);
% $$$ 
% $$$ 
% $$$ % $$$ %debugging code
% $$$ % $$$ possibleAssignmentsBPG = zeros(size(MASTER,1),size(MASTER,2));
% $$$ % $$$ for i = 1:size(MASTER,1)
% $$$ % $$$   possibleAssignmentsBPG(i,i) = 1;
% $$$ % $$$ end
% $$$ % $$$ 
% $$$ % $$$ 
% $$$ % $$$ [assignmentAccuracy, assignments] = computeAssignmentAccuracy(peakIDs, ...
% $$$ % $$$ 						  RESNUMS, possibleAssignmentsBPG);
% $$$ % $$$ 
% $$$ % $$$ prunedAssigntable      = NVR_NOE2PROB(possibleAssignmentsBPG, NOES, IALLDISTS, NTH, ...
% $$$ % $$$ 				 ROWIN, COLIN);
% $$$ % $$$ [prunedPositionRow, prunedPositionColumn] = find(prunedAssigntable ...
% $$$ % $$$ 						 - possibleAssignmentsBPG);
% $$$ % $$$ assert (isempty(prunedPositionRow));
% $$$ % $$$ %debugging code ends
% $$$ 
% $$$ %manual running code begins. to be used manually for now.
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$   DO_INDIVIDUAL_PIECE_MANUALLY = 0; DO_MANUAL_MERGING = 0; 
% $$$   STOCHASTIC_SEARCH      = 0;
% $$$   
% $$$   %individual pieces done here
% $$$   if (DO_INDIVIDUAL_PIECE_MANUALLY)
% $$$     doIndividualPieceManually(MASTER, voter, numVoters, ASSIGNTABLE, peakIndices, ...
% $$$ 			      residueIndices, NOES, ALLDISTS, IALLDISTS, NTH, ...
% $$$ 			      peakIDs, RESNUMS, RDC1, RDC2, VECTORS, MB_ShiftScore, ...
% $$$ 			      MB_RDC_Score, numCS, numRDC, HSQCDATA, csCoefficient, rdcCoefficient, noeCoefficient,numHN_NOES);
% $$$     
% $$$     keyboard
% $$$     
% $$$   elseif (DO_MANUAL_MERGING)
% $$$     %merging done here
% $$$     
% $$$     
% $$$     doMergingManually(peakIndices, ...
% $$$ 		      residueIndices, voter,numVoters,ASSIGNTABLE,NOES,ALLDISTS,IALLDISTS,NTH,peakIDs,RESNUMS,RDC1,RDC2,VECTORS, MB_ShiftScore, ...
% $$$ 		      MB_RDC_Score, numCS, numRDC, HSQCDATA, csCoefficient, rdcCoefficient, noeCoefficient,numHN_NOES);
% $$$     keyboard
% $$$   elseif (STOCHASTIC_SEARCH)
% $$$     stochasticSearch;
% $$$   else
% $$$     [foundMASTERs, finalScores, assignmentAccuracies, totalNumAssignments]=divideAndConquer(MASTER, voter, numVoters, ASSIGNTABLE, peakIndices, ...
% $$$ 						  residueIndices, NOES, ALLDISTS, IALLDISTS, NTH, ...
% $$$ 						  peakIDs, RESNUMS, RDC1, RDC2, VECTORS, MB_ShiftScore, ...
% $$$ 						  MB_RDC_Score, numCS, numRDC, HSQCDATA,csCoefficient, ...
% $$$ 						  rdcCoefficient, ...
% $$$ 						  noeCoefficient,numHN_NOES);
% $$$     
% $$$   end
% $$$   %manual running code ends.
% $$$   
% $$$   fprintf(1, 'finished divide and conquer\n');
% $$$   overallMASTER = zeros(totalNumPeaks, totalNumResidues);
% $$$   for i = 1:totalNumAssignments
% $$$     overallMASTER = overallMASTER + foundMASTERs{i};
% $$$   end
% $$$   
% $$$   if (totalNumAssignments > 0)
% $$$     aggregateMASTER = findAggregateMASTER(overallMASTER, ...
% $$$ 					  totalNumAssignments);
% $$$   else
% $$$     aggregateMASTER = MASTER;
% $$$   end
% $$$   
% $$$   
% $$$   %computeAssignmentAccuracy(peakIDs, RESNUMS, aggregateMASTER,1);
% $$$   numCorrect = 0; numAssignedPeaks = 0;
% $$$   for peakIndex = 1:size(aggregateMASTER,1)
% $$$     residueIndex = find(aggregateMASTER(peakIndex,:));
% $$$     assert (length(residueIndex) <= 1);
% $$$     if (peakIndex == residueIndex)
% $$$       numCorrect = numCorrect + 1;
% $$$     end
% $$$     if (~isempty(residueIndex))
% $$$       numAssignedPeaks = numAssignedPeaks + 1;
% $$$     end
% $$$   end
% $$$   
% $$$   fprintf(1, 'the aggregate assignment has %d correct out of %d positions.\n',numCorrect, numAssignedPeaks);
% $$$   
% $$$   numCorrectPerIter = [numCorrectPerIter numCorrect];
% $$$   numTotalAssignmentsPerIter = [numTotalAssignmentsPerIter numAssignedPeaks];
% $$$   
% $$$   %fprintf(1, 'running A* assign...\n');
% $$$   
% $$$   %[foundMASTERs, finalScores, assignmentAccuracies, totalNumAssignments] ...
% $$$   %    = a_star_assign(MASTER, voter, numVoters, ASSIGNTABLE, peakIndices, ...
% $$$   %		    residueIndices, NOES, ALLDISTS, IALLDISTS, NTH, ...
% $$$   %		    peakIDs, RESNUMS, RDC1, RDC2, VECTORS);
% $$$   
% $$$   %fprintf(1, 'aggregated the assignments.\n');
% $$$   MASTER = aggregateMASTER;
% $$$ 
% $$$   if (secondRun == 1)
% $$$   else
% $$$     S1                               = updateTen(MASTER,RDC1,VECTORS);
% $$$     S2                               = updateTen(MASTER,RDC2,VECTORS);
% $$$   end
% $$$     
% $$$ %  keyboard
% $$$   
% $$$   numAssignedPeaks = sum(sum(aggregateMASTER));
% $$$   if ((numAssignedPeaks == prevNumAssignedPeaks) | (sum(sum(MASTER)) ...
% $$$ 						    == totalNumPeaks))
% $$$     if (secondRun == 1)
% $$$       break;
% $$$     else
% $$$       secondRun = 1;
% $$$       MASTER    = initialMASTER;
% $$$       prevNumAssignedPeaks = sum(sum(initialMASTER));
% $$$       ASSIGNTABLE = initialASSIGNTABLE;
% $$$       ROWIN = initialROWIN;
% $$$       COLIN = initialCOLIN;
% $$$       SSCP = initialSSCP;
% $$$       SXCP = initialSXCP;
% $$$       CP = initialCP;
% $$$       TP = initialTP;
% $$$       HDE = initialHDE; 
% $$$       save alignmentTensorsAfterOneRoundOfCompleteAssignments.mat ...
% $$$ 	  S1 S2 S1_full S2_full;
% $$$       continue;
% $$$     end
% $$$   else
% $$$     prevNumAssignedPeaks = numAssignedPeaks;
% $$$   end
% $$$ 
% $$$ 
% $$$   fprintf(1, 'will try to remove assigned peaks and residues.\n');
% $$$   keyboard
% $$$   if (sum(sum(MASTER)) ~= totalNumPeaks)
% $$$     
% $$$ % $$$   [S1,flag]                        = updateTen(MASTER,RDC1,VECTORS);
% $$$ % $$$   assert (~flag);
% $$$ % $$$   [S2,flag]                        = updateTen(MASTER,RDC2,VECTORS);
% $$$ % $$$   assert (~flag);
% $$$ % $$$   
% $$$ % $$$   voter                            = initializeVoters(voter,SSCP,SXCP,CP,TP,HDE, RP1, RP2);
% $$$     
% $$$     [ASSIGNTABLE, voter, numVoters, ROWIN, COLIN,SSCP,SXCP,CP,TP,HDE] = removeAssignedPeaksAndResidues(MASTER, ASSIGNTABLE, voter, numVoters,peakIndices,residueIndices,SSCP,SXCP,CP,TP,HDE);
% $$$   
% $$$ % $$$   SSCP = voter{1}; SXCP = voter{2}; CP = voter{3}; TP = voter{4}; 
% $$$ % $$$   HDE = voter{5}; RP1= voter{6}; RP2 = voter{7};
% $$$     
% $$$   elseif (secondRun == 1)
% $$$     notCompletedAssignments  = 0;
% $$$   end
% $$$ end
% $$$ 
% $$$ figure; 
% $$$ plot(1:length(percentiles1), percentiles1*100, '*-', 1: ...
% $$$      length(percentiles1),percentiles2*100, 'ro-',1: ...
% $$$      length(percentiles1),numCorrectPerIter, 'kx-', 1:length(percentiles1), numTotalAssignmentsPerIter, 'g^-',1:length(percentiles1),100*numCorrectPerIter./numTotalAssignmentsPerIter,'cv-');
% $$$ grid
% $$$ legend('medium 1 MTC', 'medium 2 MTC', 'numCorrectAssignments', ...
% $$$        'numTotalAssignments','assignment accuracy','Location','BestOutside');
% $$$ 
% $$$ % $$$ [overallMASTER,foundAnAssignment, ...
% $$$ % $$$  stopGeneratingAssignments, minScore, totalNumAssignments, foundMASTERs, finalScores, assignmentAccuracies]  = assign(totalNumAssignments, finalScores, foundMASTERs, foundAnAssignment, minScore, MASTER, scoreSoFar, voter, ...
% $$$ % $$$ 						  numVoters, ASSIGNTABLE, peakIndices, residueIndices, NOES, ...
% $$$ % $$$ 						  ALLDISTS, IALLDISTS, NTH, overallMASTER, ...
% $$$ % $$$ 						  peakIDs, RESNUMS, ...
% $$$ % $$$ 						  RDC1, RDC2, VECTORS, assignmentAccuracies);
% $$$ 
% $$$ %analyzeOverallAssignments(overallMASTER, peakIDs, RESNUMS);
% $$$ 
% $$$ %keyboard
% $$$ 
% $$$ %function [SSCP,SXCP,CP,TP,HDE, RP1, RP2, ASSIGNTABLE] = removeAssignedRowsAndColumns(SSCP,SXCP,CP,TP,HDE, RP1, RP2, ...
% $$$ %						  ASSIGNTABLE, ...
% $$$ %						  relPeakIndices, relResidueIndices)
% $$$ 
