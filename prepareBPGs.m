function [RP1_A,RP2_A,HDE_A,TP_A,CP_A,SXCP_A,SSCP_A] = prepareBPGs (RDC1,VECTORS,S1,RDC2,S2,HDEXCHANGE,HBOND,peakIDs,HSHIFTS,NSHIFTS,TYPES,SSTRUCT,NOES,IALLDISTS,ALLDISTS, SHIFTS_Filename,SHIFTX_Filename,MASTER)

%given alignment tensors S1 and S2, as well as the experimental
%data, computes the bipartite graphs corresponding to each of the
%voters, to be used in the score computation.


minval      = 10e-40;

NTH         = 5.5;

ASSIGNTABLE = MASTER*0+1;

ROWIN       = 1:size(ASSIGNTABLE,1);

COLIN       = 1:size(ASSIGNTABLE,2);

RP1_A       = NVR_RDC2PROB  (ASSIGNTABLE,RDC1,VECTORS,S1)+minval;RP1_A =ren(RP1_A);

RP2_A       = NVR_RDC2PROB  (ASSIGNTABLE,RDC2,VECTORS,S2)+minval;RP2_A =ren(RP2_A);

HDE_A       = HD_HD2PROB    (ASSIGNTABLE ,HDEXCHANGE,HBOND)+minval;HDE_A=ren(HDE_A);

TP_A        = NVR_TOCSY2PROB(peakIDs,HSHIFTS,NSHIFTS,TYPES,SSTRUCT,NOES,IALLDISTS,NTH,ROWIN,COLIN)+minval;TP_A =ren(TP_A);

CP_A        = HD_CS2PROB    (ASSIGNTABLE ,HSHIFTS,NSHIFTS,TYPES,SSTRUCT,NOES,ALLDISTS,NTH,ROWIN,COLIN)+minval;CP_A =ren(CP_A);

SXCP_A      = HD_SHIFTX2PROB(ASSIGNTABLE,HSHIFTS,NSHIFTS,TYPES,SSTRUCT, ...
			     NOES,ALLDISTS,NTH,ROWIN,COLIN, ...
			     SHIFTX_Filename)+minval;SXCP_A =ren(SXCP_A);

SSCP_A      = HD_SHIFTS2PROB(ASSIGNTABLE,HSHIFTS,NSHIFTS,TYPES,SSTRUCT, ...
			     NOES,ALLDISTS,NTH,ROWIN,COLIN, SHIFTS_Filename, ...
			     SHIFTX_Filename)+minval;SSCP_A =ren(SSCP_A);
