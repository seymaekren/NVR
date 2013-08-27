dbstop if error
dbstop if warning

  if (~exist('HSQCDATA'))
  
    useOrigData     = 1;
    
    %[NOES, RDC1, RDC2, HDEXCHANGE, peakIDs, HSHIFTS, ...
    % NSHIFTS, HN_SHIFTS, RDCs] = readNMR_Data(useOrigData);
    [HSQCDATA, NOES] = loadUbqData(useOrigData);
    
    if (useOrigData)
      
      modeIndex = 7; modelIndex = 6;
      
      SHIFTX_Filename   = sprintf('InputFiles/SHIFTX/MySHIFTX.%d.model%d',modeIndex, modelIndex);
      
      SHIFTS_Filename   = sprintf('InputFiles/SHIFTS/MySHIFTS.%d.model%d',modeIndex, modelIndex);
      
%      modelDataFilename = sprintf('ModelDataGeneration/ModelDataFiles/Mode%d.coords%d',  modeIndex, modelIndex);

%      [VECTORS,TYPES,RESNUMS,SSTRUCT, HBOND, ALLDISTS,IALLDISTS] = loadmodeldata(modelDataFilename);
 
      [dontUseVECTORS,TYPES,RESNUMS,SSTRUCT, HBOND, ALLDISTS,IALLDISTS, ignoredHSQCDATA] = loaddata('InputFiles/myinput.m');

         CH_RDCS = load('InputFiles/C-H_medium1.m');
% $$$ % $$$ % $$$       %      CH_RDCS  = load('NjpData/N-H_medium2.m');
         NH_RDCS = load('InputFiles/N-H_medium1.m');
         VECTORS_CH = load('InputFiles/C-H_vectors.m');
         VECTORS_NH = load('InputFiles/N-H_vectors.m');

% $$$ 
% $$$         CH_RDCS = load('NjpData/C-H_medium1.m');
% $$$ % $$$ % $$$       %      CH_RDCS  = load('NjpData/N-H_medium2.m');
% $$$         NH_RDCS = load('NjpData/N-H_medium1.m');
% $$$         VECTORS_CH = load('NjpData/C-H_vectors.m');
% $$$         VECTORS_NH = load('NjpData/N-H_vectors.m');
      
% $$$       CH_RDCS = load('hSRI/C-H_medium1.m');
% $$$       %      CH_RDCS  = load('NjpData/N-H_medium2.m');
% $$$       NH_RDCS = load('hSRI/N-H_medium1.m');
% $$$       VECTORS_CH = load('hSRI/C-H_vectors.m');
% $$$       VECTORS_NH = load('hSRI/N-H_vectors.m');

% $$$         CH_RDCS = load('MZ_ubq/C-H_medium1.m');
% $$$ % $$$ % $$$       %      CH_RDCS  = load('NjpData/N-H_medium2.m');
% $$$         NH_RDCS = load('MZ_ubq/N-H_medium1.m');
% $$$         VECTORS_CH = load('MZ_ubq/C-H_vectors.m');
% $$$         VECTORS_NH = load('MZ_ubq/N-H_vectors.m');
% $$$       
      
      
% $$$       CH_RDCS = load('MZ_ubi/C-H_medium1.m');
% $$$ % $$$ % $$$       %      CH_RDCS  = load('NjpData/N-H_medium2.m');
% $$$       NH_RDCS = load('MZ_ubi/N-H_medium1.m');
% $$$       VECTORS_CH = load('MZ_ubi/C-H_vectors.m');
% $$$       VECTORS_NH = load('MZ_ubi/N-H_vectors.m');
      
% $$$        CH_RDCS = load('Ubi_debug/C-H_medium1.m');
% $$$ % $$$ % $$$ % $$$ % $$$ % $$$       %      CH_RDCS  = load('NjpData/N-H_medium2.m');
% $$$        NH_RDCS = load('Ubi_debug/N-H_medium1.m');
% $$$        VECTORS_CH = load('Ubi_debug/C-H_vectors.m');
% $$$        VECTORS_NH = load('Ubi_debug/N-H_vectors.m');
      
% $$$          CH_RDCS = load('Ubi_debug2/C-H_medium1.m');
% $$$ % $$$ % $$$ % $$$ % $$$ % $$$       %      CH_RDCS  = load('NjpData/N-H_medium2.m');
% $$$          NH_RDCS = load('Ubi_debug2/N-H_medium1.m');
% $$$          VECTORS_CH = load('Ubi_debug2/C-H_vectors.m');
% $$$          VECTORS_NH = load('Ubi_debug2/N-H_vectors.m');

      CH_RDCS = CH_RDCS(:,2);
      NH_RDCS = NH_RDCS(:,2);
      VECTORS_CH = VECTORS_CH(:,2:4);
      VECTORS_NH = VECTORS_NH(:,2:4);
      fprintf(1, 'loaded 1ubi vectors and RDCs.\n');
      

    else
      %use distrib data
      
      assert (0);
      
      [VECTORS,TYPES,RESNUMS,SSTRUCT, HBOND, ALLDISTS,IALLDISTS] = ...
	  loadmodeldata('modeldata.txt');
      
    end
    
    useMBM_EM = 1;
    
    %betterHD(HSQCDATA(:,2:3),HSQCDATA(:,4:5),HSQCDATA(:,6),HSQCDATA(:,1),NOES,VECTORS, ...
    % 	 TYPES,RESNUMS,SSTRUCT, HBOND, ALLDISTS,IALLDISTS, ...
    % 	 SHIFTS_Filename, SHIFTX_Filename, useMBM_EM);
    
    
    peaks      = HSQCDATA(:,2:3);
 %   rdcs       = HSQCDATA(:,4:5);
    rdcs       = [NH_RDCS, CH_RDCS];
    HDEXCHANGE = HSQCDATA(:,6) ;
    peakIDs    = HSQCDATA(:,1) ;
  end

betterNVR_withCH


%[oldNVR_SCORE, assignmentAccuracy, assignments, NVR_SCORE]= NVR_withCH(peaks,NH_RDCS, CH_RDCS,HDEXCHANGE, peakIDs, NOES, VECTORS_NH, VECTORS_CH,TYPES, ...
%						  RESNUMS,SSTRUCT, HBOND, ALLDISTS, SHIFTS_Filename, SHIFTX_Filename);

%testOneMediumAlignmentTensor(NH_RDCS, CH_RDCS,VECTORS_NH, VECTORS_CH);


%[oldNVR_SCORE, assignmentAccuracy, assignments, NVR_SCORE]= NVR(peaks,rdcs,HDEXCHANGE, peakIDs, NOES, VECTORS_NH,TYPES, ...
%						  RESNUMS,SSTRUCT, HBOND, ALLDISTS, SHIFTS_Filename, SHIFTX_Filename);


%fprintf(1, 'the assignment accuracy is %f\n', assignmentAccuracy);


% $$$ [assignmentAccuracy,NVR_SCORE,oldNVR_SCORE]= betterNVR(HSQCDATA(:,2:3),HSQCDATA(:,4:5),HSQCDATA(:,6),HSQCDATA(:, 1),NOES,VECTORS, ...
% $$$ 						  TYPES,RESNUMS,SSTRUCT, HBOND, ALLDISTS,SHIFTS_Filename, SHIFTX_Filename);



