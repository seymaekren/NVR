fprintf(1, 'clearing environment variables...\n');
clear all;

b_running1CMZ                     = 0;
b_runningMBP                      = 0;

cd InputFiles/

system('./checkFiles.sh')
cd ..
    
useCH_RDCs                    = 0;

refineWithRDCs                  = 0;
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
    
    
    peaks      = HSQCDATA(:,2:3);
    rdcs       = HSQCDATA(:,4:5);
    HDEXCHANGE = HSQCDATA(:,6) ;
    peakIDs    = HSQCDATA(:,1) ;
  end

  if (b_printDistances)
    printDistances(ALLDISTS);
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


