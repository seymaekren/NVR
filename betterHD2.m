function betterHD2(peaks, rdcs, HDEXCHANGE, NH_RDCS, CH_RDCS, ...
    peakIDs, NOES, VECTORS_NH, VECTORS_CH, VECTORS, TYPES, ...
    RESNUMS,SSTRUCT, HBOND, ALLDISTS, IALLDISTS, SHIFTS_Filename, ...
    SHIFTX_Filename, useCH_RDCs,useHD_Routines, useTOCSY, ...
    refineWithRDCs, truncateProbabilities, ...
    b_printOptimizationFiles, b_printIndividualScores, useFourRDCsPerResidue, ...
    b_runningMBP, b_runningEIN, b_running1FQB, b_runningPoln, ...
    b_printIndividualScoringMatrix, inputModeIndex, inputModelIndex)

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
    
    elseif ((refineWithRDCs == 1) && (b_runningMBP == 1))
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
    
    else %changed
        [ROWIN, COLIN, ASSIGNTABLE, MASTER, HDE, TP, CP, SXCP, SSCP, RDC1, RDC2, NTH, ...
            RP1, RP2, S1, S2] = initialize(peaks,rdcs,...
            HDEXCHANGE, ...
            peakIDs, NOES, ...
            VECTORS,TYPES, ...
            RESNUMS,SSTRUCT, HBOND, ...
            ALLDISTS,IALLDISTS, ...
            SHIFTS_Filename, ...
            SHIFTX_Filename, ...
            useCH_RDCs,useHD_Routines, useTOCSY, truncateProbabilities, ...
            b_runningMBP, b_runningEIN, b_running1FQB,b_runningPoln); 
    end
 
    b_reduceSizeOfTheProblem = 0;

    if (b_reduceSizeOfTheProblem)  
        reduceSizeOfTheProblem();
    end
    
    if (b_printOptimizationFiles) || (b_printIndividualScores) 
        
        myMASTER                  = MASTER*0;
        numPeaks                  = size(myMASTER,1);
        
        if (refineWithRDCs)
        
            %changed
            fprintf('burdaaaa');
            %myMASTER_WithExtraColumns =     load (assignmentMatrixFilename);
            myMASTER_WithExtraColumns = eye(70,75);
            fprintf(1, 'the matrix read has %d rows and %d columns.\n', ...
                size(myMASTER_WithExtraColumns,1), size(myMASTER_WithExtraColumns,2));
            myMASTER                  =     myMASTER_WithExtraColumns(:,1:size(myMASTER_WithExtraColumns,2)-3);
            fprintf(1, 'the derived matrix has %d columns.\n', size(myMASTER,2));
          
                         
            %changed-end
            
            %             myMASTER_WithExtraColumns =     load (assignmentMatrixFilename);
            %             fprintf(1, 'reading %s...\n', assignmentMatrixFilename);
            %             fprintf(1, 'the matrix read has %d rows and %d columns.\n', ...
            %                 size(myMASTER_WithExtraColumns,1), size(myMASTER_WithExtraColumns,2));
            %             myMASTER                  =     myMASTER_WithExtraColumns(:,1:size(myMASTER_WithExtraColumns,2)-3);
            %             fprintf(1, 'the derived matrix has %d columns.\n', size(myMASTER,2));

            numCorrectAssignmentsInTheFileRead = 0;
            
            for peakIndex = 1:numPeaks
                if (myMASTER(peakIndex,peakIndex) == 1)
                    numCorrectAssignmentsInTheFileRead = ...
                        numCorrectAssignmentsInTheFileRead + 1;
                end
            end
            fprintf(1, 'assignment accuracy in the inputted file = %f\n',numCorrectAssignmentsInTheFileRead/size(myMASTER,1));

            if (useCH_RDCs == 1)
                S1                               = updateTen_CH(myMASTER, NH_RDCS,CH_RDCS, VECTORS_NH, VECTORS_CH);
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
                    S1                               = updateTen(myMASTER,CCa_RDCS,VECTORS_CCa);
                    %S2                               = updateTen(myMASTER,NH_RDCS,VECTORS_NH)
                    S3                               = updateTen(myMASTER,NC_RDCS,VECTORS_NC);
                    
                    RP1                              = NVR_RDC2PROB(ASSIGNTABLE,CCa_RDCS,VECTORS_CCa,S1,ROWIN,COLIN);
                    %RP2                              = NVR_RDC2PROB(ASSIGNTABLE,NH_RDCS,VECTORS_NH,S2,ROWIN,COLIN);
                    RP3                              = NVR_RDC2PROB(ASSIGNTABLE,NC_RDCS,VECTORS_NC,S3,ROWIN,COLIN);
                else
                    S1                               = updateTen(myMASTER,RDC1,VECTORS);
                    S2                               = updateTen(myMASTER,RDC2,VECTORS);
                    RP1                              = NVR_RDC2PROB(ASSIGNTABLE,RDC1,VECTORS,S1,ROWIN,COLIN);
                    RP2                              = NVR_RDC2PROB(ASSIGNTABLE,RDC2,VECTORS,S2,ROWIN,COLIN);
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
            
            elseif ((b_runningMBP == 1) && (refineWithRDCs == 1))
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
                elseif  (b_printCS)
                    scoringMatrixFilename = 'individualScoringMatrix_CS.txt';
                    printOptimizationFiles_IndScoringMatrix(CP, NOES, ALLDISTS, ...
                        COLIN, NTH, scoringMatrixFilename, SXCP, SSCP);
                elseif  (b_printRDC)
                    assert (refineWithRDCs == 1);
                    scoringMatrixFilename = 'individualScoringMatrix_RDC.txt';
                    printOptimizationFiles_IndScoringMatrix(RP1, NOES, ALLDISTS, ...
                        COLIN, NTH, scoringMatrixFilename, RP2);
                end
            else
                
                %changed
                printOptimizationFiles(CP, SXCP, SSCP, TP, HDE, refineWithRDCs, ...
                    RP1, RP2, NOES, ALLDISTS, COLIN, NTH, ...
                    useBayesianMatrix);
            end
        
        elseif (b_printIndividualScores)
            printIndividualScores(CP, SXCP, SSCP, TP, HDE, refineWithRDCs, ...
                RP1, RP2);
            filterScores
        end
    
    end
end
