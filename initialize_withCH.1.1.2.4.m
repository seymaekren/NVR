function [ROWIN, COLIN, ASSIGNTABLE, MASTER, HDE, CP, SXCP, SSCP, NTH, RP1, RP2] = initialize_withCH(peaks,HDEXCHANGE, peakIDs, NOES, VECTORS_NH, VECTORS_CH,TYPES, RESNUMS,SSTRUCT, ...
						  HBOND, ALLDISTS, ...
						  IALLDISTS, ...
						  SHIFTS_Filename, ...
						  SHIFTX_Filename, ...
						  NH_RDCS, CH_RDCS);

HSHIFTS = peaks(:,1);
NSHIFTS = peaks(:,2);
%dontUseRDC1 = rdcs(:,1);
%dontUseRDC2 = rdcs(:,2);
RDC1    = NH_RDCS;
RDC2    = CH_RDCS;

%compute a tollerance for NOE distances
NTH=4.8;
mu=mean(mean(ALLDISTS));
if(mu-12.9>0)
   NTH=NTH+(mu-12.9); NTH=min(NTH,8);
end
fprintf(1, 'the original NTH is %f\n',NTH);
NTH = 8.60; 		
%NTH = 6.43; 
%NTH =  8.48;
fprintf(1, 'set NTH to %f Angstroms.\n',NTH);
keyboard
ASSIGNTABLE = ones(length(HSHIFTS),size(VECTORS_NH,1))/size(VECTORS_NH,1);
%these keep track of which peaks and residues are represented in the current
%matricies
ROWIN       = 1:size(ASSIGNTABLE,1);
COLIN       = 1:size(ASSIGNTABLE,2);
%This is the master assignment table
MASTER      = ASSIGNTABLE*0;

%HDE = NVR_HD2PROB(ASSIGNTABLE,HDEXCHANGE,HBOND);
HDE = ones(size(ASSIGNTABLE,1),size(ASSIGNTABLE,2));
fprintf(1, 'bypassing using HD data.\n');
%fprintf(1,'check out HDE matrix. type return to cont.\n');
keyboard

ASSIGNTABLE = and(ASSIGNTABLE,HDE);
 nlast = sum(sum(ASSIGNTABLE));
 for(i=1:100)
    [ASSIGNTABLE] = lockdown(MASTER,ASSIGNTABLE,NOES,ALLDISTS, NTH,ROWIN,COLIN);
    if(sum(sum(ASSIGNTABLE )) == nlast)
       break;
    end
    nlast = sum(sum(ASSIGNTABLE ));
 end

[CP] = NVR_CS2PROB(ASSIGNTABLE,HSHIFTS,NSHIFTS,TYPES,SSTRUCT,NOES,ALLDISTS,NTH,ROWIN,COLIN);SXCP=CP;SSCP = CP;RP1 = CP;RP2 = CP;
save CP.ubi.debug.mat CP ASSIGNTABLE HSHIFTS NSHIFTS TYPES SSTRUCT NOES ...
    ALLDISTS NTH ROWIN COLIN
fprintf(1, 'saved CP.ubi.debug.mat\n');
keyboard

ASSIGNTABLE = and(ASSIGNTABLE,CP);
fprintf(1, 'CP loaded.\n');
%keyboard
 for(lv = [round(size(CP,2)*.6)])
    V = vote(CP,CP,CP,CP,CP,ASSIGNTABLE,NOES,ALLDISTS,NTH,ROWIN,COLIN,2,lv);
    V2 = vote(CP',CP',CP',CP',CP',ASSIGNTABLE',NOES,ALLDISTS,NTH,ROWIN,COLIN,2,lv);
    V = or(V,V2');
    ASSIGNTABLE = ASSIGNTABLE.*V(1:size(ASSIGNTABLE,1),1:size(ASSIGNTABLE,2));
    
    nlast = sum(sum(ASSIGNTABLE));
    for(i=1:100)
       [ASSIGNTABLE] = lockdown(MASTER,ASSIGNTABLE,NOES,ALLDISTS, NTH,ROWIN,COLIN);
       if(sum(sum(ASSIGNTABLE)) == nlast)
          break;
       end
       nlast = sum(sum(ASSIGNTABLE));
    end
 end


%prune via the program shiftx
[SXCP]      = NVR_SHIFTX2PROB(ASSIGNTABLE,HSHIFTS,NSHIFTS,TYPES,SSTRUCT, ...
			 NOES,ALLDISTS,NTH,ROWIN,COLIN, SHIFTX_Filename);
ASSIGNTABLE = and(ASSIGNTABLE,SXCP);
fprintf(1, 'SXCP loaded.\n');
%keyboard
 for(lv = [round(size(CP,2)*.5)])
    %reduce to a fixed percentage
    V = vote(SXCP,SXCP,SXCP,SXCP,SXCP,ASSIGNTABLE,NOES,ALLDISTS,NTH,ROWIN,COLIN,2,lv);
    V2 = vote(SXCP',SXCP',SXCP',SXCP',SXCP',ASSIGNTABLE',NOES,ALLDISTS,NTH,ROWIN,COLIN,2,lv);
    V = or(V,V2');
    ASSIGNTABLE = ASSIGNTABLE.*V(1:size(ASSIGNTABLE,1),1:size(ASSIGNTABLE,2));
    nlast = sum(sum(ASSIGNTABLE));
    for(i=1:100)
       [ASSIGNTABLE] = lockdown(MASTER,ASSIGNTABLE,NOES,ALLDISTS, NTH,ROWIN,COLIN);
       if(sum(sum(ASSIGNTABLE)) == nlast)
          break;
       end
       nlast = sum(sum(ASSIGNTABLE));
    end
 end


%prune via the program shifts
[SSCP] = NVR_SHIFTS2PROB(ASSIGNTABLE,HSHIFTS,NSHIFTS,TYPES,SSTRUCT, ...
			 NOES,ALLDISTS,NTH,ROWIN,COLIN, SHIFTS_Filename, SHIFTX_Filename);


ASSIGNTABLE = and(ASSIGNTABLE,SSCP);

fprintf(1, 'SSCP loaded.\n');
%keyboard





function V = vote(CP,SXCP,SSCP,RP1,RP2,ASSIGNTABLE,NOES,ALLDISTS,NTH,ROWIN,COLIN,atype,top)
%singles
[A]=getASS(CP,CP, ASSIGNTABLE,0,atype,top);
[B]=getASS(SXCP,SXCP, ASSIGNTABLE,0,atype,top);
[C]=getASS(SSCP,SSCP, ASSIGNTABLE,0,atype,top);
[D]=getASS(RP1,RP1, ASSIGNTABLE,0,atype,top);
[E]=getASS(RP2,RP2, ASSIGNTABLE,0,atype,top);
%doubles
[F]=getASS(SXCP.*CP,SXCP.*CP, ASSIGNTABLE,0,atype,top);
[G]=getASS(SXCP.*RP1,SXCP.*RP1, ASSIGNTABLE,0,atype,top);
[H]=getASS(SXCP.*RP2,SXCP.*RP2, ASSIGNTABLE,0,atype,top);
[I]=getASS(SXCP.*SSCP,SXCP.*SSCP, ASSIGNTABLE,0,atype,top);
[J]=getASS(SSCP.*CP,SSCP.*CP, ASSIGNTABLE,0,atype,top);
[K]=getASS(SSCP.*RP1,SSCP.*RP1, ASSIGNTABLE,0,atype,top);
[L]=getASS(SSCP.*RP2,SSCP.*RP2, ASSIGNTABLE,0,atype,top);
[M]=getASS(CP.*RP1,CP.*RP1, ASSIGNTABLE,0,atype,top);
[N]=getASS(CP.*RP2,CP.*RP2, ASSIGNTABLE,0,atype,top);
[O]=getASS(RP1.*RP2,RP1.*RP2, ASSIGNTABLE,0,atype,top);
%triples
[P]=getASS(CP.*SXCP.*SSCP,CP.*SXCP.*SSCP, ASSIGNTABLE,0,atype,top);
[Q]=getASS(CP.*SXCP.*RP1,CP.*SXCP.*RP1, ASSIGNTABLE,0,atype,top);
[R]=getASS(CP.*SXCP.*RP2,CP.*SXCP.*RP2, ASSIGNTABLE,0,atype,top);
[S]=getASS(CP.*SSCP.*RP1,CP.*SSCP.*RP1, ASSIGNTABLE,0,atype,top);
[T]=getASS(CP.*SSCP.*RP2,CP.*SSCP.*RP2, ASSIGNTABLE,0,atype,top);
[U]=getASS(SXCP.*SSCP.*RP1,SXCP.*SSCP.*RP1, ASSIGNTABLE,0,atype,top);
[V]=getASS(SXCP.*SSCP.*RP2,SXCP.*SSCP.*RP2, ASSIGNTABLE,0,atype,top);
[W]=getASS(SXCP.*RP1.*RP2,SXCP.*RP1.*RP2, ASSIGNTABLE,0,atype,top);
[X]=getASS(SSCP.*RP1.*RP2,SSCP.*RP2.*RP1, ASSIGNTABLE,0,atype,top);
[Y]=getASS(CP.*RP1.*RP2,CP.*RP2.*RP1, ASSIGNTABLE,0,atype,top);
%quads
[Z]=getASS(CP.*RP1.*RP2.*SXCP,CP.*RP2.*RP1.*SXCP, ASSIGNTABLE,0,atype,top);
[AA]=getASS(CP.*RP1.*RP2.*SSCP,CP.*RP2.*RP1.*SSCP, ASSIGNTABLE,0,atype,top);
[BB]=getASS(SXCP.*RP1.*RP2.*SSCP,SXCP.*RP2.*RP1.*SSCP, ASSIGNTABLE,0,atype,top);
[CC]=getASS(SXCP.*CP.*RP2.*SSCP,SXCP.*RP2.*SSCP.*CP, ASSIGNTABLE,0,atype,top);
[DD]=getASS(SXCP.*CP.*RP1.*SSCP,SXCP.*RP1.*SSCP.*CP, ASSIGNTABLE,0,atype,top);
[EE]=getASS(RP1.*SXCP.*CP.*RP2.*SSCP,RP1.*SXCP.*RP2.*SSCP.*CP, ASSIGNTABLE,0,atype,top);
V = A+B+C+D+E+F+G+H+I+J+K+L+M+N+O+P+Q+R+S+T+U+V+W+X+Y+Z+AA+BB+CC+DD+EE;



