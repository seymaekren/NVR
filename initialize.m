function [ROWIN, COLIN, ASSIGNTABLE, MASTER, HDE, TP, CP, SXCP, SSCP, ...
	  RDC1, RDC2, NTH, RP1, RP2, S1, S2] ...
    = initialize(peaks, rdcs, HDEXCHANGE, peakIDs, NOES, ...
		 VECTORS,TYPES, RESNUMS,SSTRUCT, ...
		 HBOND, ALLDISTS, ...
		 IALLDISTS, ...
		 SHIFTS_Filename, SHIFTX_Filename, useCH_RDCs, ...
		 useHD_Routines, useTOCSY, truncateProbabilities, ...
		 b_runningMBP, b_runningEIN, b_running1FQB,b_runningPoln)

%truncateProbabilities is only taken into account if NVR routines
%are used.

if (useTOCSY)
  fprintf(1, 'using TOCSY data file...\n');
else
  fprintf(1, 'not using TOCSY data file ...\n');
end

if (useHD_Routines)
  fprintf(1, 'using HD routines, with their thresholding and NOE applying...\n');
else
  fprintf(1, 'using NVR routines, with their smaller range of deviation...\n');
end
HSHIFTS = peaks(:,1);
NSHIFTS = peaks(:,2);
CASHIFTS = peaks(:,3);
RDC1    = rdcs(:,1);
RDC2    = rdcs(:,2);

%compute a tollerance for NOE distances
NTH=4.8;
mu=mean(mean(ALLDISTS));
if(mu-12.9>0)
   NTH=NTH+(mu-12.9); 
   %   NTH=min(NTH,8);
   
end


NTH = NTH + 2.34; %this is due to 1G6J max. violation distance.

NTH=min(NTH,9.33);


fprintf(1, 'NTH = %f\n', NTH);

ASSIGNTABLE = ones(length(HSHIFTS),size(VECTORS,1))/size(VECTORS,1);
OASSIGNTABLE=ASSIGNTABLE;
%these keep track of which peaks and residues are represented in the current
%matricies
ROWIN=1:size(ASSIGNTABLE,1);
COLIN=1:size(ASSIGNTABLE,2);
%This is the master assignment table
MASTER=ASSIGNTABLE*0;

fprintf('computing assignments...\n');
if (useHD_Routines)
  HDE         = myHD2PROB(ASSIGNTABLE,HDEXCHANGE,HBOND);
else
  HDE= NVR_HD2PROB(ASSIGNTABLE,HDEXCHANGE,HBOND, truncateProbabilities);
end
ASSIGNTABLE = and(ASSIGNTABLE,HDE);

foundZeroEntryOnDiagonal = 0;
for i=1:min(size(ASSIGNTABLE,1),size(ASSIGNTABLE,2))
   if (ASSIGNTABLE(i,i) == 0)
     fprintf(1, 'after AND with HDE, i = %d ASSIGNTABLE(i,i) = 0\n',i);
     foundZeroEntryOnDiagonal = 1;
%     keyboard
   end
end

TP =   ASSIGNTABLE;
%fprintf(1, 'warning. TOCSY data not used. this is good for ff2, hSRI.\n');
%keyboard
if (useTOCSY)
  TP = NVR_TOCSY2PROB(peakIDs,HSHIFTS,NSHIFTS,TYPES,SSTRUCT,NOES, ...
		      IALLDISTS,ROWIN,COLIN);
  
end
ASSIGNTABLE = and(ASSIGNTABLE,TP );



for i=1:min(size(ASSIGNTABLE,1),size(ASSIGNTABLE,2))
   if (ASSIGNTABLE(i,i) == 0)
     fprintf(1, 'after AND with TP, i = %d ASSIGNTABLE(i,i) = 0\n', ...
	     i);
     foundZeroEntryOnDiagonal = 1;
   end
end

if (useHD_Routines)
  [CP] = HD_CS2PROB(OASSIGNTABLE,HSHIFTS,NSHIFTS,TYPES,SSTRUCT, ...
		    NOES,IALLDISTS,NTH,ROWIN,COLIN);
else
  [CP] = NVR_CS2PROB(ASSIGNTABLE,HSHIFTS,NSHIFTS,CASHIFTS,TYPES,SSTRUCT, ...
		     NOES,ALLDISTS,NTH,ROWIN,COLIN, truncateProbabilities, ...
		     b_runningMBP, b_runningEIN, b_running1FQB);
end
ASSIGNTABLE = and(ASSIGNTABLE,CP);

for i=1:min(size(ASSIGNTABLE,1),size(ASSIGNTABLE,2))
   if (ASSIGNTABLE(i,i) == 0)
     fprintf(1, 'after AND with CP, i = %d ASSIGNTABLE(i,i) = 0\n',i);
     foundZeroEntryOnDiagonal = 1;
   end
end


%prune via the program shiftx

if (useHD_Routines)
  [SXCP] = HD_SHIFTX2PROB(OASSIGNTABLE,HSHIFTS,NSHIFTS,CASHIFTS,TYPES,SSTRUCT, ...
			  NOES,IALLDISTS,NTH,ROWIN,COLIN, SHIFTX_Filename);
else
  [SXCP]      = NVR_SHIFTX2PROB(ASSIGNTABLE,HSHIFTS,NSHIFTS,CASHIFTS,TYPES,SSTRUCT, ...
				NOES,ALLDISTS,NTH,ROWIN,COLIN, ...
				SHIFTX_Filename, truncateProbabilities, ...
				b_runningMBP, b_runningEIN, ...
				b_runningPoln, b_running1FQB);
end
ASSIGNTABLE = and(ASSIGNTABLE,SXCP);   

for i=1:min(size(ASSIGNTABLE,1),size(ASSIGNTABLE,2))
   if (ASSIGNTABLE(i,i) == 0)
     fprintf(1, 'After AND with SXCP, i = %d ASSIGNTABLE(i,i) = 0\n',i);
     foundZeroEntryOnDiagonal = 1;
     %     keyboard
   end
end


%prune via the program shifts
% if (useHD_Routines)
%   [SSCP] = HD_SHIFTS2PROB(OASSIGNTABLE,HSHIFTS,NSHIFTS,TYPES,CASHIFTS,SSTRUCT, ...
% 			  NOES,IALLDISTS,NTH,ROWIN,COLIN, SHIFTS_Filename, ...
% 			  SHIFTX_Filename);
% else
%   [SSCP] = NVR_SHIFTS2PROB(ASSIGNTABLE,HSHIFTS,NSHIFTS,CASHIFTS,TYPES,SSTRUCT, ...
% 			   NOES,ALLDISTS,NTH,ROWIN,COLIN, SHIFTS_Filename, ...
% 			   SHIFTX_Filename, truncateProbabilities, ...
% 			   b_runningMBP, b_runningEIN, b_running1FQB);
% end
% ASSIGNTABLE = and(ASSIGNTABLE,SSCP);
% 
% for i=1:min(size(ASSIGNTABLE,1),size(ASSIGNTABLE,2))
%    if (ASSIGNTABLE(i,i) == 0)
%      fprintf(1, 'After AND with SSCP, i = %d ASSIGNTABLE(i,i) = 0\n',i);
%      foundZeroEntryOnDiagonal = 1;
%    end
% end
% 
% if (foundZeroEntryOnDiagonal)
%   fprintf(1, 'found zero entries on diagonal of the matrices. proceed?\n');
%   keyboard
% end

SSCP=ASSIGNTABLE;
%SXCP=ASSIGNTABLE;
%SSCP=SSCP.*ASSIGNTABLE;
SXCP=SXCP.*ASSIGNTABLE;
CP=CP.*ASSIGNTABLE;
TP=TP.*ASSIGNTABLE;
HDE=HDE.*ASSIGNTABLE;

S1 = []; S2 = []; RP1 = []; RP2 = []; 

%if (useCH_RDCs)
%   fprintf(1, 'initialize assigned the alignment tensor');
%   fprintf(1, ' and score matrices to 0.\n');
%   fprintf(1, 'later can compute the correct matrices and tensors.\n');
%else
%  S1 = updateTen(MASTER,RDC1,VECTORS);
%  S2 = updateTen(MASTER,RDC2,VECTORS);
%  RP1 = NVR_RDC2PROB(ASSIGNTABLE,RDC1,VECTORS,S1, ROWIN, COLIN);
%  RP2 = NVR_RDC2PROB(ASSIGNTABLE,RDC2,VECTORS,S2, ROWIN, COLIN);
%end
%MASTER = origMASTER;

% $$$ marsShiftScore = computeMarsShiftScore(ASSIGNTABLE, HSHIFTS, NSHIFTS, ...
% $$$ 				       SHIFTS_Filename);
% $$$ 
% $$$ [MB_ShiftScore,numCS]  = computeMB_ShiftScore(ASSIGNTABLE, HSHIFTS, NSHIFTS, ...
% $$$ 				       SHIFTS_Filename);

