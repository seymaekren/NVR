fprintf(1, 'clearing environment variables...\n');
clear all;

b_running1CMZ                     = 0;
b_runningMBP                      = 0;

cd InputFiles/

% $$$ !./cleanSymbolicLinks.sh
% $$$ 
% $$$ if (b_running1CMZ)
% $$$   !./symbolicLinksFor1CMZ.sh
% $$$ else
% $$$   !./symbolicLinksForCAM2.sh
% $$$ end

    system("./checkFiles.sh")
cd ..
    
fprintf(1, 'warning. file setting is commented out\n');
keyboard

useCH_RDCs                    = 0;

refineWithRDCs                = 1;

%assignmentMatrixFilename      = '/home2/apaydin/Workdir/OptimizationFiles/1D3Z/1UBI/WithNVR_TOCSY/WithNTH=8.80/NOEsFromMZ_ubq/TruncatingWithLargerThresholds/WithoutRDC/1UBI.txt';
%assignmentMatrixFilename      = '/home2/apaydin/Workdir/OptimizationFiles/1D3Z/1UBI/WithNVR_TOCSY/WithNTH=9.33/WithRDC_FirstRound/1UBI.txt';
%assignmentMatrixFilename      = '/home2/apaydin/Workdir/OptimizationFiles/1D3Z/1UBQ/WithNVR/WithoutRDCs/WithTOCSY/WithMZ_NOE/TruncatingWithLargerThreshold/NTH=8.82/1UBQ.txt';
%assignmentMatrixFilename      = '/home2/apaydin/Workdir/OptimizationFiles/1D3Z/1UD7/WithNVR/WithoutRDCs/WithTOCSY/MZ_NOEs/NTH=8.71/TruncatingWithLargerThreshold/1UD7.txt';
%assignmentMatrixFilename      = '/home2/apaydin/Workdir/OptimizationFiles/1D3Z/1G6J/WithNVR/WithoutRDCs/WithTOCSY/TruncatingWithLargerThreshold/WithMZ_NOE/NTH=8.65/1G6J.txt';
%assignmentMatrixFilename      = '/home2/apaydin/Workdir/OptimizationFiles/hSRI/WithNVR_AndTOCSY/WithHD-Exchange/WithNTH=9.33/WithoutTruncating/WithoutTOCSY_Thresholding/WithoutRDCs/hSRI.txt';
%assignmentMatrixFilename      = '/home2/apaydin/Workdir/Weighting/hSRI/WithRDCsFirstRound/hSRI.txt';
%assignmentMatrixFilename       = '/home2/apaydin/Workdir/OptimizationFiles/ff2/WithNVR_AndTOCSY/WithNTH=9.33/WithNewTruncation/WithoutRDC/ff2.txt';
%assignmentMatrixFilename      = '/home2/apaydin/Workdir/OptimizationFiles/3GB1/1GB1/WithoutRDCs/WithNTH7.14/TruncationWithLargerThreshold/1GB1.txt';
%assignmentMatrixFilename      = '/home2/apaydin/Workdir/OptimizationFiles/3GB1/2GB1/WithoutRDC/WithNTH=7.14/TruncationWithLargerThreshold/2GB1.txt';
%assignmentMatrixFilename      = '/home2/apaydin/Workdir/OptimizationFiles/3GB1/1PGB/WithoutRDC/WithTOCSY/WithNewTruncation/WithNTH=7.27/1PGB.txt';
%assignmentMatrixFilename      = '/home2/apaydin/Workdir/OptimizationFiles/1E8L/193L/WithTOCSY/WithNTH=9.33/TruncatingWithNewBounds/WithoutRDC/193L.txt';
%assignmentMatrixFilename      = '/home2/apaydin/Workdir/OptimizationFiles/1E8L/193L/WithTOCSY/WithNTH=9.33/TruncatingWithNewBounds/WithRDC_SecondRound/193L.txt';
%assignmentMatrixFilename      = '/home2/apaydin/Workdir/OptimizationFiles/1E8L/1AKI/WithTOCSY/TruncatingWithNewBounds/NTH=9.33/WithoutRDCs/1AKI.txt';
%assignmentMatrixFilename      = '/home2/apaydin/Workdir/OptimizationFiles/1E8L/1AZF/WithoutRDC/1AZF.txt';
%assignmentMatrixFilename      = '/home2/apaydin/Workdir/Weighting/1E8L/1AZF/WithoutRDCs/1AZF.txt';
%assignmentMatrixFilename      = '/home2/apaydin/Workdir/OptimizationFiles/1E8L/1BGI/WithoutRDC/1BGI.txt';
%assignmentMatrixFilename      = '/home2/apaydin/Workdir/OptimizationFiles/1E8L/1H87/WithoutRDCs/1H87.txt';
%assignmentMatrixFilename      = '/home2/apaydin/Workdir/OptimizationFiles/1E8L/1LSC/WithoutRDCs/1LSC.txt';
%assignmentMatrixFilename       = '/home2/apaydin/Workdir/OptimizationFiles/1E8L/1LSE/WithoutRDCs/1LSE.txt';
%assignmentMatrixFilename       ='/home2/apaydin/Workdir/OptimizationFiles/1E8L/2LYZ/WithoutRDCs/2LYZ.txt';
%assignmentMatrixFilename       = '/home2/apaydin/Workdir/OptimizationFiles/1E8L/3LYZ/WithoutRDCs/3LYZ.txt';
%assignmentMatrixFilename       = '/home2/apaydin/Workdir/OptimizationFiles/1E8L/5LYZ/WithoutRDCs/5LYZ.txt';
%assignmentMatrixFilename       = '/home2/apaydin/Workdir/OptimizationFiles/1E8L/6LYZ/WithoutRDCs/6LYZ.txt';
%assignmentMatrixFilename       = '/home2/apaydin/Workdir/OptimizationFiles/1D3Z/Control/1ESR/1ESR.txt';
%assignmentMatrixFilename       = '/home2/apaydin/Workdir/OptimizationFiles/1D3Z/1AAR/WithRDC/WithHD-exchange/With4RDCsPerResidue/1AAR.txt';
%assignmentMatrixFilename       = '/home2/apaydin/Workdir/OptimizationFiles/Poln/WithTOCSY/WithHD-Exchange/POLN.txt';
%assignmentMatrixFilename        = '/home2/apaydin/Workdir/OptimizationFiles/GB1/WithTOCSY/WithHD-exchange/WithRDC/1GB1.txt';
%assignmentMatrixFilename        = '/home2/apaydin/Workdir/OptimizationFiles/GB1/WithTOCSY/WithHD-exchange/WithRDC/1GB1.txt';
%assignmentMatrixFilename        = '/home2/apaydin/Workdir/OptimizationFiles/1D3Z/1UBI/WithNVR_TOCSY/WithNTH=8.80/NOEsFromMZ_ubq/TruncatingWithLargerThresholds/WithRDC_FirstRound/WithCH_and_NH_RDCs/1UBI.txt';
%assignmentMatrixFilename        = '/home2/apaydin/Workdir/OptimizationFiles/1CMZ/1CMZ/WithoutRDC/1CMZ.txt';
%assignmentMatrixFilename        = '/home2/apaydin/Workdir/MBP_100PercentCorrect.txt';
assignmentMatrixFilename        = '/home2/apaydin/Workdir/1LYZ_100PercentCorrect.txt';
b_printOptimizationFiles        = 1;
useHD_Routines                  = 0; 
useTOCSY                        = 1;
runOriginalNVR                  = 0;
runOriginalHD                   = 0;
b_printIndividualScores         = 0; %this is outputted to obtain the
                                   %individual scores for weighting
useFourRDCsPerResidue           = 0;
b_printNumAvailableRDC_Values   = 0;
b_printRP1                      = 0;
b_printRP2                      = 0;
b_printHDE                      = 0;

b_printIndividualScoringMatrix  = 0; %this is outputted to test the
				     %contribution of each data
				     %source to the assignments
b_printCP                       = 0;
b_printSSCP                     = 0;
b_printSXCP                     = 0;
b_printTP                       = 0;
b_printCS                       = 0;
b_printRDC                      = 0;
b_printDistances                = 1;

truncateProbabilities           = 1; %this parameter only taken into
				       %account if using NVR routines, not
				       %HD routines, and in NVR_CS2PROB,
				       %NVR_SHIFTS2PROB, NVR_SHIFTX2PROB.



				       
				       
if (truncateProbabilities)
  fprintf(1, 'truncating probabilities...\n');
else
  fprintf(1, 'not truncating probabilities...\n');
end

if (~useTOCSY)
  fprintf(1, 'useTOCSY is 0. Please enter return\n');
  keyboard
end

if (useFourRDCsPerResidue)
  assert (useCH_RDCs == 1);
  assert (refineWithRDCs == 1);
  fprintf(1, 'using 4 RDCs/residue\n');
end
      
if (b_printIndividualScoringMatrix)
  fprintf(1, 'will print individual scoring matrix.\n');
  assert (b_printOptimizationFiles == 1);
end

if (b_printIndividualScores)
  fprintf(1, 'will print individual scores in libsvm format.\n');
  assert (b_printOptimizationFiles);
end


dbstop if error
dbstop if warning

  if (~exist('HSQCDATA'))
  
    useOrigData      = 1;
    [HSQCDATA, NOES] = readNMR_Data2(useOrigData);
    
%    hsqc_fid = fopen('hsqcData_ff2.txt', 'w');
%    for hsqc_index = 1:size(HSQCDATA,1)
%      fprintf(hsqc_fid, '%f %f\n', HSQCDATA(hsqc_index,2), HSQCDATA(hsqc_index,3));
%    end
%    fclose(hsqc_fid);
%    fprintf(1, 'printed HSQC data into hsqcData_ff2.txt\n');
%    keyboard
    
    if (useOrigData)
      
      modeIndex = 7; modelIndex = 6;
      
      SHIFTX_Filename   = sprintf('InputFiles/SHIFTX/MySHIFTX.%d.model%d',modeIndex, modelIndex);
      
      SHIFTS_Filename   = sprintf('InputFiles/SHIFTS/MySHIFTS.%d.model%d',modeIndex, modelIndex);
      
      if (useCH_RDCs)

	 CH_RDCS    = load('InputFiles/C-H_medium1.m');
         NH_RDCS    = load('InputFiles/N-H_medium1.m');
         VECTORS_CH = load('InputFiles/C-H_vectors.m');
         VECTORS_NH = load('InputFiles/N-H_vectors.m');

	 CH_RDCS    = CH_RDCS(:,2);
	 NH_RDCS    = NH_RDCS(:,2);
	 VECTORS_CH = VECTORS_CH(:,2:4);
	 VECTORS_NH = VECTORS_NH(:,2:4);

	 fprintf(1, 'total CH RDCs:%d\n', length(CH_RDCS));
	 fprintf(1, 'total NH RDCs:%d\n', length(NH_RDCS));
	 
	 if (b_printNumAvailableRDC_Values)
	 
	   printNumAvailableRDC_Values(CH_RDCS, NH_RDCS);
	   
	 end
	 
	 fprintf(1, 'loaded CH vectors and RDCs.\n');
	 [doNotUseVECTORS,TYPES,RESNUMS,SSTRUCT, HBOND, ALLDISTS,IALLDISTS, ...
	  ignoredHSQCDATA] = loaddata('InputFiles/myinput.m');

       elseif ((refineWithRDCs == 1) & (b_runningMBP == 1))

         CCa_RDCS    = load('InputFiles/C-Ca_medium1.m');
         NC_RDCS     = load('InputFiles/N-C_medium1.m');
         NH_RDCS     = load('InputFiles/N-H_medium1.m');
         VECTORS_CCa = load('InputFiles/C-Ca_vectors.m');
         VECTORS_NC  = load('InputFiles/N-C_vectors.m');
         VECTORS_NH  = load('InputFiles/N-H_vectors.m');

         CCa_RDCS    = CCa_RDCS(:,2);
         NC_RDCS     = NC_RDCS(:,2);
         NH_RDCS     = NH_RDCS(:,2);
         VECTORS_CCa = VECTORS_CCa(:,2:4);
         VECTORS_NH  = VECTORS_NH(:,2:4);
         VECTORS_NC  = VECTORS_NC(:,2:4);


         fprintf(1, 'total CCa RDCs:%d\n', length(CCa_RDCS));
         fprintf(1, 'total NH RDCs:%d\n', length(NH_RDCS));
         fprintf(1, 'total NC RDCs:%d\n', length(NC_RDCS));

         if (b_printNumAvailableRDC_Values)

           printNumAvailableRDC_Values3(CCa_RDCS, NH_RDCS, NC_RDCS);

         end

         fprintf(1, 'loaded vectors and RDCs.\n');
         [doNotUseVECTORS,TYPES,RESNUMS,SSTRUCT, HBOND, ALLDISTS,IALLDISTS, ...
	  ...
          ignoredHSQCDATA] = loaddata('InputFiles/myinput.m');

      else

	 [VECTORS,TYPES,RESNUMS,SSTRUCT, HBOND, ALLDISTS,IALLDISTS, ...
	  ignoredHSQCDATA] = loaddata('InputFiles/myinput.m');
	 
	 if (b_running1CMZ)
	   VECTORS_NH = load('InputFiles/N-H_vectors.m');
	   VECTORS    = VECTORS_NH(:,2:4);
	   fprintf(1, 'running 1CMZ and overwriting 1DK8 vectors with 1CMZ vectors.\n');
	 end
      
      end

    else
      %use distrib data
      
      assert (0);
      
      [VECTORS,TYPES,RESNUMS,SSTRUCT, HBOND, ALLDISTS,IALLDISTS] = ...
	  loadmodeldata('modeldata.txt');
      
    end
    
    %betterHD(HSQCDATA(:,2:3),HSQCDATA(:,4:5),HSQCDATA(:,6),HSQCDATA(:,1),NOES,VECTORS, ...
    % 	 TYPES,RESNUMS,SSTRUCT, HBOND, ALLDISTS,IALLDISTS, ...
    % 	 SHIFTS_Filename, SHIFTX_Filename, useMBM_EM);
    
    
    peaks      = HSQCDATA(:,2:3);
    rdcs       = HSQCDATA(:,4:5);
    HDEXCHANGE = HSQCDATA(:,6) ;
    peakIDs    = HSQCDATA(:,1) ;
%    fprintf(1, 'peak ids are read here.\n');
%    keyboard
  end

  if (b_printDistances)
    printDistances(ALLDISTS);
    fprintf(1, 'continue? Enter return if yes, dbquit if not\n');
    keyboard
  end

  if ((~runOriginalNVR) & (~runOriginalHD))
    betterHD
  elseif (runOriginalNVR)
    
    if (useCH_RDCs)
      [oldNVR_SCORE, assignmentAccuracy, assignments, NVR_SCORE]= NVR_withCH(peaks,NH_RDCS, CH_RDCS,HDEXCHANGE, peakIDs, NOES, VECTORS_NH, VECTORS_CH,TYPES, ...
	     RESNUMS,SSTRUCT, HBOND, ALLDISTS, IALLDISTS, SHIFTS_Filename, ...
						  SHIFTX_Filename, truncateProbabilities);
    else
      [oldNVR_SCORE, assignmentAccuracy, assignments, NVR_SCORE]= NVR(peaks,rdcs,HDEXCHANGE, peakIDs, NOES, VECTORS,TYPES, ...
						  RESNUMS,SSTRUCT, ...
						  HBOND, ALLDISTS, ...
						  IALLDISTS, ...
						  SHIFTS_Filename, ...
						  SHIFTX_Filename, ...
						  truncateProbabilities);
    end
    
  elseif (runOriginalHD)
    [HD_SCORE,assignmentAccuracy, assignments, weightedHD_Score] = HD(peaks,rdcs,HDEXCHANGE, ...
						  peakIDs, NOES, ...
						  VECTORS,TYPES, ...
						  RESNUMS,SSTRUCT, HBOND, ALLDISTS, ...
						  IALLDISTS, SHIFTS_Filename, SHIFTX_Filename)
    
  end
  %keyboard
%fprintf(1, 'the assignment accuracy is %f\n', assignmentAccuracy);


% $$$ [assignmentAccuracy,NVR_SCORE,oldNVR_SCORE]= betterNVR(HSQCDATA(:,2:3),HSQCDATA(:,4:5),HSQCDATA(:,6),HSQCDATA(:, 1),NOES,VECTORS, ...
% $$$ 						  TYPES,RESNUMS,SSTRUCT, HBOND, ALLDISTS,SHIFTS_Filename, SHIFTX_Filename);



