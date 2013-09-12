function runBetter (inputModeIndex, inputModelIndex)


fprintf(1, 'clearing environment variables...\n');
%clear all;

b_running1CMZ                     = 0;
b_runningMBP                      = 0;

%this should be changed 

cd InputFiles/
    system('./checkFiles.sh')
cd ../

fprintf(1, 'warning. file setting is commented out\n');

useCH_RDCs                    = 0;

refineWithRDCs                = 0;

b_printOptimizationFiles        = 1;
useHD_Routines                  = 0; 
useTOCSY                        = 1;
runOriginalNVR                  = 0;
runOriginalHD                   = 0;
b_printIndividualScores         = 0; %this is outputted to obtain the individual scores for weighting
useFourRDCsPerResidue           = 0;
b_printNumAvailableRDC_Values   = 0;
b_printRP1                      = 0;
b_printRP2                      = 0;
b_printHDE                      = 0;

b_printIndividualScoringMatrix  = 0; %this is outputted to test the contribution of each data
                                        %source to the assignments
b_printCP                       = 0;
b_printSSCP                     = 0;
b_printSXCP                     = 0;
b_printTP                       = 0;
b_printCS                       = 0;
b_printRDC                      = 0;
b_printDistances                = 1;

truncateProbabilities = 1; %this parameter only taken into account if using NVR routines, not
                            %HD routines, and in NVR_CS2PROB, NVR_SHIFTS2PROB, NVR_SHIFTX2PROB.

%simple printouts
simplePrintOuts(truncateProbabilities, useTOCSY, useFourRDCsPerResidue,...
    useCH_RDCs, refineWithRDCs, b_printIndividualScoringMatrix, b_printOptimizationFiles, b_printIndividualScores);

dbstop if error
dbstop if warning

[NOES, NH_RDCS, CH_RDCS, VECTORS_NH, VECTORS_CH, VECTORS, TYPES,RESNUMS,SSTRUCT, HBOND, ALLDISTS, IALLDISTS,...
    peaks, rdcs, HDEXCHANGE, peakIDs, SHIFTS_Filename,SHIFTX_Filename] ...
    = loadHSQCdata(useCH_RDCs, b_printNumAvailableRDC_Values, refineWithRDCs, b_runningMBP, b_running1CMZ, inputModeIndex, inputModelIndex);


printDistances(b_printDistances, ALLDISTS);

runNVRorHD(runOriginalNVR, runOriginalHD, useCH_RDCs, ...
    peaks, rdcs, NH_RDCS, CH_RDCS,HDEXCHANGE, peakIDs, NOES, VECTORS_NH, VECTORS_CH, VECTORS, TYPES, ...
    RESNUMS,SSTRUCT, HBOND, ALLDISTS, IALLDISTS, SHIFTS_Filename, SHIFTX_Filename, ... 
    truncateProbabilities, useHD_Routines, useTOCSY, refineWithRDCs, ...
    b_printOptimizationFiles, b_printIndividualScores, useFourRDCsPerResidue, ...
    b_runningMBP, b_printIndividualScoringMatrix, inputModeIndex, inputModelIndex);
  
        
        
        
%%%%%%%%%%%%%%%%%%%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%  

function simplePrintOuts(truncateProbabilities, useTOCSY, useFourRDCsPerResidue,...
    useCH_RDCs, refineWithRDCs, b_printIndividualScoringMatrix, b_printOptimizationFiles, b_printIndividualScores)

if (truncateProbabilities)
    fprintf(1, 'truncating probabilities...\n');
else
    fprintf(1, 'not truncating probabilities...\n');
end

if (~useTOCSY)
    fprintf(1, 'useTOCSY is 0. Please enter return\n');
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [NOES, NH_RDCS, CH_RDCS, VECTORS_NH, VECTORS_CH, VECTORS, TYPES, RESNUMS, SSTRUCT, HBOND, ALLDISTS, IALLDISTS,...
    peaks, rdcs, HDEXCHANGE, peakIDs, SHIFTS_Filename,SHIFTX_Filename] ...
    = loadHSQCdata(useCH_RDCs, b_printNumAvailableRDC_Values, refineWithRDCs, b_runningMBP, b_running1CMZ, inputModeIndex, inputModelIndex)

useOrigData      = 1;

%HSQCDATA****************************************
[HSQCDATA, NOES] = readNMR_Data2(useOrigData);%changed

if (useOrigData)
    
    if (inputModeIndex == 78)
        filename       = sprintf('InputFiles/1EF1/ModelDataGeneration/ModelDataFiles/Mode.7.8.coords%d',  inputModelIndex);
    else
        filename       = sprintf('InputFiles/1EF1/ModelDataGeneration/ModelDataFiles/Mode%d.coords%d',  inputModeIndex, inputModelIndex);
    end  
    
    SHIFTX_Filename   = sprintf('InputFiles/1EF1/SHIFTX/MySHIFTX.%d.model%d',inputModeIndex, inputModelIndex);
    SHIFTS_Filename   = sprintf('InputFiles/1EF1/SHIFTS/MySHIFTS.%d.model%d',inputModeIndex, inputModelIndex);




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
        [doNotUseVECTORS,TYPES,RESNUMS,SSTRUCT, HBOND, ALLDISTS,IALLDISTS] = loadmodeldata(filename);  
        

    elseif ((refineWithRDCs == 1) && (b_runningMBP == 1)) 

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
        [doNotUseVECTORS,TYPES,RESNUMS,SSTRUCT, HBOND, ALLDISTS,IALLDISTS] = loadmodeldata(filename);  

    else 

        [VECTORS,TYPES,RESNUMS,SSTRUCT, HBOND, ALLDISTS, IALLDISTS] = loadmodeldata(filename);
        
        NH_RDCS = [];
        CH_RDCS = [];
        VECTORS_NH = []; 
        VECTORS_CH = [];
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function printDistances(b_printDistances, ALLDISTS)

if (b_printDistances)
    printDistances(ALLDISTS);
    fprintf(1, 'continue? Enter return if yes, dbquit if not\n');
end


function runNVRorHD(runOriginalNVR, runOriginalHD, useCH_RDCs, ...
    peaks, rdcs, NH_RDCS, CH_RDCS,HDEXCHANGE, peakIDs, NOES, VECTORS_NH, VECTORS_CH, VECTORS, TYPES, ...
    RESNUMS,SSTRUCT, HBOND, ALLDISTS, IALLDISTS, SHIFTS_Filename, ...
    SHIFTX_Filename, truncateProbabilities, useHD_Routines, useTOCSY, ...
    refineWithRDCs, b_printOptimizationFiles, b_printIndividualScores, useFourRDCsPerResidue, ...
    b_runningMBP, b_printIndividualScoringMatrix, inputModeIndex, inputModelIndex)

if ((~runOriginalNVR) && (~runOriginalHD)) 
    betterHD(peaks, rdcs, HDEXCHANGE, NH_RDCS, CH_RDCS, ...
        peakIDs, NOES, VECTORS_NH, VECTORS_CH, VECTORS, TYPES, ...
        RESNUMS,SSTRUCT, HBOND, ALLDISTS, IALLDISTS, SHIFTS_Filename, ...
        SHIFTX_Filename, useCH_RDCs,useHD_Routines, useTOCSY, ...
        refineWithRDCs, truncateProbabilities, ...
        b_printOptimizationFiles, b_printIndividualScores, useFourRDCsPerResidue, ...
        b_runningMBP, b_printIndividualScoringMatrix, inputModeIndex, inputModelIndex); %changed

elseif (runOriginalNVR)
    if (useCH_RDCs)
        [oldNVR_SCORE, assignmentAccuracy, assignments, NVR_SCORE]= NVR_withCH(peaks,NH_RDCS, CH_RDCS,HDEXCHANGE, peakIDs, NOES, VECTORS_NH, VECTORS_CH, TYPES, ...
            RESNUMS,SSTRUCT, HBOND, ALLDISTS, IALLDISTS, SHIFTS_Filename, ...
            SHIFTX_Filename, truncateProbabilities);
    else
        [oldNVR_SCORE, assignmentAccuracy, assignments, NVR_SCORE]= NVR(peaks,rdcs,HDEXCHANGE, peakIDs, NOES, VECTORS, TYPES, ...
            RESNUMS,SSTRUCT, HBOND, ALLDISTS, IALLDISTS, SHIFTS_Filename, ...
            SHIFTX_Filename, truncateProbabilities);
    end
    
elseif (runOriginalHD)
    [HD_SCORE,assignmentAccuracy, assignments, weightedHD_Score] = HD(peaks,rdcs,HDEXCHANGE, ...
        peakIDs, NOES, VECTORS, TYPES, RESNUMS,SSTRUCT, HBOND, ALLDISTS, ...
        IALLDISTS, SHIFTS_Filename, SHIFTX_Filename);
end

