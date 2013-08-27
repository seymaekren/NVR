function [MASTER, ASSIGNTABLE, CP, SXCP, SSCP,  HDE, ROWIN, COLIN] = ...
    AssignFirstFiveResidues_withCH(MASTER, ASSIGNTABLE, NOES, IALLDISTS, ...
			    NTH, ROWIN, COLIN, ALLDISTS, CP, SXCP, ...
			    SSCP, RP1, RP2, HDE, S1, NH_RDCS, ...
			    CH_RDCS, VECTORS_NH, VECTORS_CH)


if (0)
  origSSCP        = SSCP;
  origMASTER      = MASTER;
  origASSIGNTABLE = ASSIGNTABLE;
  origROWIN       = ROWIN;
  origCOLIN       = COLIN;
  origCP          = CP;
  origSSCP        = SSCP;
  origSXCP        = SXCP;
  origHDE         = HDE;
  origRP1         = RP1;
  origRP2         = RP2;
  
  
  
  
  
  
  for(lv = [round(size(CP,2)*.35)])
    %reduce to a fixed percentage
    V = vote(SSCP,SSCP,SSCP,SSCP,SSCP,ASSIGNTABLE,NOES,ALLDISTS,NTH,ROWIN,COLIN,2,lv);
    V2 = vote(SSCP',SSCP',SSCP',SSCP',SSCP',ASSIGNTABLE',NOES,ALLDISTS,NTH,ROWIN,COLIN,2,lv);
    V = or(V,V2');
    ASSIGNTABLE = ASSIGNTABLE.*V(1:size(ASSIGNTABLE,1),1:size(ASSIGNTABLE,2));
    
    
    nlast = sum(sum(ASSIGNTABLE));
    for(i=1:100)
      [MASTER,ASSIGNTABLE]=updateMASTER(MASTER,ASSIGNTABLE,ROWIN,COLIN);[ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN);
      [ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,ROWIN,COLIN]=reduceAll(ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,ROWIN,COLIN);[ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN);
      [ASSIGNTABLE] = lockdown(MASTER,ASSIGNTABLE,NOES,ALLDISTS, NTH,ROWIN,COLIN);
      [MASTER,ASSIGNTABLE]=updateMASTER(MASTER,ASSIGNTABLE,ROWIN,COLIN);[ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN);
      [ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,ROWIN,COLIN]=reduceAll(ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,ROWIN,COLIN);[ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN);
      if(sum(sum(ASSIGNTABLE)) == nlast)
	break;
      end
      
      
      nlast = sum(sum(ASSIGNTABLE));
    end
    
  end
  
  
  
  
  
  
  
  
  [ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN);
  [MASTER,ASSIGNTABLE]=updateMASTER(MASTER,ASSIGNTABLE,ROWIN,COLIN);[ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN);
  [ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,ROWIN,COLIN]=reduceAll(ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,ROWIN,COLIN);[ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN);
  [ASSIGNTABLE] = lockdown(MASTER,ASSIGNTABLE,NOES,ALLDISTS, NTH,ROWIN,COLIN);
  [MASTER,ASSIGNTABLE]=updateMASTER(MASTER,ASSIGNTABLE,ROWIN,COLIN);[ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN);
  [ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,ROWIN,COLIN]=reduceAll(ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,ROWIN,COLIN);[ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN);
  
  fprintf(1, 'assigned %d residues.\n', sum(sum(MASTER)));



if(sum(sum(MASTER))>=5)
   
  S1 = updateTen_CH(MASTER,NH_RDCS,CH_RDCS, VECTORS_NH, VECTORS_CH);
   RP1 = NVR_RDC2PROB(ASSIGNTABLE,NH_RDCS,VECTORS_NH,S1, ROWIN, COLIN);
   RP2 = NVR_RDC2PROB_CH(ASSIGNTABLE,CH_RDCS,VECTORS_CH,S1, ROWIN, COLIN);
else
   S1 = ones(3);S2=S1;
   RP1 = CP*0+1;
   for(i=1:size(RP1,1))
      RP1(i,:)=RP1(i,:)/sum(RP1(i,:));
   end
   RP2 = RP1;
end


%ok, now reduce to a linear number of edges
for(lv = [30,28,26,24,22,20,18,16])
   V = vote(CP,SXCP,SSCP,RP1,RP2,ASSIGNTABLE,NOES,ALLDISTS,NTH,ROWIN,COLIN,2,lv);
   V2 = vote(CP',SXCP',SSCP',RP1',RP2',ASSIGNTABLE',NOES,ALLDISTS,NTH,ROWIN,COLIN,2,lv);
   V = or(V,V2');
   ASSIGNTABLE = ASSIGNTABLE.*V(1:size(ASSIGNTABLE,1),1:size(ASSIGNTABLE,2));
   
   nlast = sum(sum(ASSIGNTABLE));
   for(i=1:100)
      [MASTER,ASSIGNTABLE]=updateMASTER(MASTER,ASSIGNTABLE,ROWIN,COLIN);[ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN);
      [ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,ROWIN,COLIN]=reduceAll(ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,ROWIN,COLIN);[ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN);
      [ASSIGNTABLE] = lockdown(MASTER,ASSIGNTABLE,NOES,ALLDISTS, NTH,ROWIN,COLIN);
      [MASTER,ASSIGNTABLE]=updateMASTER(MASTER,ASSIGNTABLE,ROWIN,COLIN);[ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN);
      [ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,ROWIN,COLIN]=reduceAll(ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,ROWIN,COLIN);[ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN);
      if(sum(sum(ASSIGNTABLE)) == nlast)
         break;
      end
      nlast = sum(sum(ASSIGNTABLE));
   end
end

%Phase 1: make the first few assignments
last = sum(sum(ASSIGNTABLE));
while(sum(sum(MASTER))<8)
   if(size(MASTER,1)-sum(sum(MASTER))>0)   
      V = vote3(CP,SXCP,SSCP,ASSIGNTABLE,NOES,ALLDISTS,NTH,ROWIN,COLIN,1,1);
      nV=CP.*SXCP.*SSCP;nV(1:size(ASSIGNTABLE,1),1:size(ASSIGNTABLE,2))=nV(1:size(ASSIGNTABLE,1),1:size(ASSIGNTABLE,2)).*ASSIGNTABLE;
      nV=thresh(nV,.004);
      
      V(1:size(ASSIGNTABLE,1),1:size(ASSIGNTABLE,2))=V(1:size(ASSIGNTABLE,1),1:size(ASSIGNTABLE,2)).*ASSIGNTABLE;
      [A]=getASS(V,V, ASSIGNTABLE,0,1);A(1:size(ASSIGNTABLE,1),1:size(ASSIGNTABLE,2))=A(1:size(ASSIGNTABLE,1),1:size(ASSIGNTABLE,2)).*ASSIGNTABLE;
      V=thresh(V,7);
      
      X = and(A(1:size(ASSIGNTABLE,1),1:size(ASSIGNTABLE,2)),and(V(1:size(ASSIGNTABLE,1),1:size(ASSIGNTABLE,2)),nV));
      
      for(i=1:size(X,1))
         x = find(X(i,:));
         if(length(x)==1)
            ASSIGNTABLE(i,:)=0;
            ASSIGNTABLE(i,x)=1;
         end
      end
   end
   
   ilast = sum(sum(ASSIGNTABLE));
   for(i=[.3 .3 .3 .3 .25 .25 .25 .25 .2 .2 .2 .2 .1 .1 .1 .1])
     [MASTER,ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,ROWIN,COLIN]= ...
	  combine(SXCP.*CP,SSCP,MASTER,ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,ROWIN,COLIN,NOES,ALLDISTS,NTH,S1,NH_RDCS,CH_RDCS,VECTORS_NH,VECTORS_CH, i);
      [MASTER,ASSIGNTABLE]=updateMASTER(MASTER,ASSIGNTABLE,ROWIN,COLIN);[ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN);
      [ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,ROWIN,COLIN]=reduceAll(ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,ROWIN,COLIN);[ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN);
      [ASSIGNTABLE] = lockdown(MASTER,ASSIGNTABLE,NOES,ALLDISTS, NTH,ROWIN,COLIN);
      [MASTER,ASSIGNTABLE]=updateMASTER(MASTER,ASSIGNTABLE,ROWIN,COLIN);[ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN);
      [ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,ROWIN,COLIN]=reduceAll(ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,ROWIN,COLIN);[ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN);
      
      
      
      S1 = updateTen_CH(MASTER,NH_RDCS,CH_RDCS, VECTORS_NH, VECTORS_CH);
      RP1 = NVR_RDC2PROB(ASSIGNTABLE,NH_RDCS,VECTORS_NH,S1, ROWIN, COLIN);
      RP2 = NVR_RDC2PROB_CH(ASSIGNTABLE,CH_RDCS,VECTORS_CH,S1, ROWIN, COLIN);
      
      
      if(size(MASTER,1)-sum(sum(MASTER))<=1)
         break;
      end
      if(sum(sum(ASSIGNTABLE))==ilast)
         break;
      end
      ilast = sum(sum(ASSIGNTABLE));
   end
   
   if(sum(sum(ASSIGNTABLE))==last)
      break;
   end
   last = sum(sum(ASSIGNTABLE));
   if(sum(sum(MASTER))>9)
      break;
   end
end

%assign any that are unambiguous via chem shifts
ilast = sum(sum(ASSIGNTABLE));
for(i=1:20)
  [MASTER,ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,ROWIN,COLIN]= ...
      combine(SXCP.*CP,SSCP,MASTER,ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,ROWIN,COLIN,NOES,ALLDISTS,NTH,S1,NH_RDCS,CH_RDCS,VECTORS_NH,VECTORS_CH, .9);
   [MASTER,ASSIGNTABLE]=updateMASTER(MASTER,ASSIGNTABLE,ROWIN,COLIN);[ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN);
   [ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,ROWIN,COLIN]=reduceAll(ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,ROWIN,COLIN);[ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN);
   [ASSIGNTABLE] = lockdown(MASTER,ASSIGNTABLE,NOES,ALLDISTS, NTH,ROWIN,COLIN);
   [MASTER,ASSIGNTABLE]=updateMASTER(MASTER,ASSIGNTABLE,ROWIN,COLIN);[ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN);
   [ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,ROWIN,COLIN]=reduceAll(ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,ROWIN,COLIN);[ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN);
   
         
      
   
   
   S1 = updateTen_CH(MASTER,NH_RDCS,CH_RDCS, VECTORS_NH, VECTORS_CH);
   RP1 = NVR_RDC2PROB(ASSIGNTABLE,NH_RDCS,VECTORS_NH,S1, ROWIN, COLIN);
   RP2 = NVR_RDC2PROB_CH(ASSIGNTABLE,CH_RDCS,VECTORS_CH,S1, ROWIN, COLIN);
   
   
   if(size(MASTER,1)-sum(sum(MASTER))<=1)
      break;
   end
   if(sum(sum(ASSIGNTABLE))==ilast)
      break;
   end
   ilast = sum(sum(ASSIGNTABLE));
end




S1 = updateTen_CH(MASTER,NH_RDCS,CH_RDCS, VECTORS_NH, VECTORS_CH);
RP1 = NVR_RDC2PROB(ASSIGNTABLE,NH_RDCS,VECTORS_NH,S1, ROWIN, COLIN);
RP2 = NVR_RDC2PROB_CH(ASSIGNTABLE,CH_RDCS,VECTORS_CH,S1, ROWIN, COLIN);




for(q=[.9999,.999,.99,.9])
  [MASTER,ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,ROWIN,COLIN]= ...
       combine(RP1,RP2,MASTER,ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,ROWIN,COLIN,NOES,ALLDISTS,NTH,S1,NH_RDCS,CH_RDCS,VECTORS_NH,VECTORS_CH, q);[ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN);
   [MASTER,ASSIGNTABLE]=updateMASTER(MASTER,ASSIGNTABLE,ROWIN,COLIN);[ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN);
   [ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,ROWIN,COLIN]=reduceAll(ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,ROWIN,COLIN);[ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN);
   [ASSIGNTABLE] = lockdown(MASTER,ASSIGNTABLE,NOES,ALLDISTS, NTH,ROWIN,COLIN);
   [MASTER,ASSIGNTABLE]=updateMASTER(MASTER,ASSIGNTABLE,ROWIN,COLIN);[ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN);
   [ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,ROWIN,COLIN]=reduceAll(ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,ROWIN,COLIN);[ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN);
   
   
   S1 = updateTen_CH(MASTER,NH_RDCS,CH_RDCS, VECTORS_NH, VECTORS_CH);
   RP1 = NVR_RDC2PROB(ASSIGNTABLE,NH_RDCS,VECTORS_NH,S1, ROWIN, COLIN);
   RP2 = NVR_RDC2PROB_CH(ASSIGNTABLE,CH_RDCS,VECTORS_CH,S1, ROWIN, COLIN);

   
   
   if(size(MASTER,1)-sum(sum(MASTER))<=1)
      break;
   end
end


%further reduce the number of edges
for(lv = [14 12 10 8 6 4])
   V = vote(CP,SXCP,SSCP,RP1,RP2,ASSIGNTABLE,NOES,ALLDISTS,NTH,ROWIN,COLIN,2,lv);
   V2 = vote(CP',SXCP',SSCP',RP1',RP2',ASSIGNTABLE',NOES,ALLDISTS,NTH,ROWIN,COLIN,2,lv);
   V = or(V,V2');
   ASSIGNTABLE = ASSIGNTABLE.*V(1:size(ASSIGNTABLE,1),1:size(ASSIGNTABLE,2));
   
   nlast = sum(sum(ASSIGNTABLE));
   for(i=1:100)
      [MASTER,ASSIGNTABLE]=updateMASTER(MASTER,ASSIGNTABLE,ROWIN,COLIN);[ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN);
      [ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,ROWIN,COLIN]=reduceAll(ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,ROWIN,COLIN);[ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN);
      [ASSIGNTABLE] = lockdown(MASTER,ASSIGNTABLE,NOES,ALLDISTS, NTH,ROWIN,COLIN);
      [MASTER,ASSIGNTABLE]=updateMASTER(MASTER,ASSIGNTABLE,ROWIN,COLIN);[ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN);
      [ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,ROWIN,COLIN]=reduceAll(ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,ROWIN,COLIN);[ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN);

      
      
      
      S1 = updateTen_CH(MASTER,NH_RDCS,CH_RDCS, VECTORS_NH, VECTORS_CH);
      RP1 = NVR_RDC2PROB(ASSIGNTABLE,NH_RDCS,VECTORS_NH,S1, ROWIN, COLIN);
      RP2 = NVR_RDC2PROB_CH(ASSIGNTABLE,CH_RDCS,VECTORS_CH,S1, ROWIN, COLIN);
      
      
      if(sum(sum(ASSIGNTABLE)) == nlast)
         break;
      end
      nlast = sum(sum(ASSIGNTABLE));
   end
end

[assignedPeakIndices, assignedResidueIndices] = find(MASTER-origMASTER);
scores                                        = zeros(length(assignedPeakIndices),1);
for peakIndex=1:length(assignedPeakIndices)
  assignedPeakIndex    = assignedPeakIndices   (peakIndex);
  assignedResidueIndex = assignedResidueIndices(peakIndex);
  scores(peakIndex)    = -log(origCP(assignedPeakIndex, ...
				assignedResidueIndex)) ;
  scores(peakIndex)    = scores(peakIndex) - log (origSSCP(assignedPeakIndex, ...
				assignedResidueIndex)) ;
  scores(peakIndex)    = scores(peakIndex)  - log (origSXCP(assignedPeakIndex, ...
						  assignedResidueIndex)) ;
   scores(peakIndex)    = scores(peakIndex)  - log(origHDE(assignedPeakIndex, ...
						  assignedResidueIndex)) ;
  
end

[sortedScores, relPeakIndices]        = sort(scores); 
MASTER                                = origMASTER;
ASSIGNTABLE                           = origASSIGNTABLE;
ROWIN                                 = origROWIN;
COLIN                                 = origCOLIN;
CP                                    = origCP;
SSCP                                  = origSSCP;
SXCP                                  = origSXCP;
HDE                                   = origHDE;
RP1                                   = origRP1;
RP2                                   = origRP2;
end
% $$$ 
%assignedPeakIndices = [1 2 3 4 5];
%assignedResidueIndices = [1 2 3 4 5];
assignedPeakIndices     = [1 2 4 5 6];
assignedResidueIndices  = [1 2 4 5 6];
%assignedPeakIndices    = [4 6 10 12 13];
%assignedResidueIndices = [4 6 10 12 13];

%assert (length(assignedPeakIndices)>=5);
%for i = 1:min(5, length(relPeakIndices))
%  peakIndex                             = assignedPeakIndices(relPeakIndices(i));
%  residueIndex                          = assignedResidueIndices(relPeakIndices(i));
for i = 1:length(assignedPeakIndices)
  peakIndex                             = assignedPeakIndices(i);
  residueIndex                          = assignedResidueIndices(i);
  ASSIGNTABLE (peakIndex,:)             = 0;
  ASSIGNTABLE (peakIndex, residueIndex) = 1;
%  fprintf(1, '%d th assignment has score %f\n',i, scores(relPeakIndices(i)));
end




[S1,saupeWrongFlag]                               = ...
    updateTen_CH(ASSIGNTABLE,NH_RDCS,CH_RDCS, VECTORS_NH, VECTORS_CH);
assert (saupeWrongFlag == 0);
%[S2,saupeWrongFlag]                               = updateTen(ASSIGNTABLE,RDC2,VECTORS);
assert (saupeWrongFlag == 0);
RP1                                               = NVR_RDC2PROB(ASSIGNTABLE,NH_RDCS,VECTORS_NH,S1,ROWIN,COLIN);
RP2                                               = NVR_RDC2PROB(ASSIGNTABLE,CH_RDCS,VECTORS_CH,S1,ROWIN,COLIN);



nlast = sum(sum(ASSIGNTABLE));
for(i=1:100)
  [MASTER,ASSIGNTABLE]=updateMASTER(MASTER,ASSIGNTABLE,ROWIN,COLIN);[ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN);
  [ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,ROWIN,COLIN]=reduceAll(ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,ROWIN,COLIN);[ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN);
% $$$   [ASSIGNTABLE] = lockdown(MASTER,ASSIGNTABLE,NOES,ALLDISTS, NTH,ROWIN,COLIN);
% $$$   [MASTER,ASSIGNTABLE]=updateMASTER(MASTER,ASSIGNTABLE,ROWIN,COLIN);[ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN);
% $$$   [ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,ROWIN,COLIN]=reduceAll(ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,ROWIN,COLIN);[ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN);
  if(sum(sum(ASSIGNTABLE)) == nlast)
    break;
  end
  nlast = sum(sum(ASSIGNTABLE));
end






function [ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,ROWIN,COLIN]=reduceAll(ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,ROWIN,COLIN);
[ASSIGNTABLE,CP,x,y]=reduce(ASSIGNTABLE,CP,ROWIN,COLIN);
[ASSIGNTABLE,SXCP,x,y]=reduce(ASSIGNTABLE,SXCP,ROWIN,COLIN);
[ASSIGNTABLE,SSCP,x,y]=reduce(ASSIGNTABLE,SSCP,ROWIN,COLIN);
[ASSIGNTABLE,ASSIGNTABLE,ROWIN,COLIN]=reduce(ASSIGNTABLE,ASSIGNTABLE,ROWIN,COLIN);
ASSIGNTABLE=and(ASSIGNTABLE,ASSIGNTABLE);



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













function V = vote3(CP,SXCP,SSCP,ASSIGNTABLE,NOES,ALLDISTS,NTH,ROWIN,COLIN,atype,top)
%singles
[A]=getASS(CP,CP, ASSIGNTABLE,0,atype,top);
[B]=getASS(SXCP,SXCP, ASSIGNTABLE,0,atype,top);
[C]=getASS(SSCP,SSCP, ASSIGNTABLE,0,atype,top);
%doubles
[D]=getASS(SXCP.*CP,SXCP.*CP, ASSIGNTABLE,0,atype,top);
[E]=getASS(SXCP.*SSCP,SXCP.*SSCP,  ASSIGNTABLE,0,atype,top);
[F]=getASS(SSCP.*CP,SSCP.*CP,  ASSIGNTABLE,0,atype,top);
%triples
[G]=getASS(CP.*SXCP.*SSCP,CP.*SXCP.*SSCP,  ASSIGNTABLE,0,atype,top);
V = A+B+C+D+E+F+G;

function  [MASTER,ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,ROWIN,COLIN]=combine(A,B,MASTER,ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,ROWIN,COLIN,NOES,ALLDISTS,NTH,S1,RDC1,RDC2,VECTORS,VECTORS_CH,THR);
THR;
if(THR>0)
   [ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN);
   [ASSIGNTABLE] = getUACMB(A,B,ASSIGNTABLE,NOES,ALLDISTS,NTH,THR, ...
			      1,ROWIN,COLIN, MASTER);[ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN);
   [MASTER,ASSIGNTABLE]=updateMASTER(MASTER,ASSIGNTABLE,ROWIN,COLIN);[ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN);
   [ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,ROWIN,COLIN]=reduceAll(ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,ROWIN,COLIN);[ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN);
   [ASSIGNTABLE] = lockdown(MASTER,ASSIGNTABLE,NOES,ALLDISTS, NTH,ROWIN,COLIN);
   [MASTER,ASSIGNTABLE]=updateMASTER(MASTER,ASSIGNTABLE,ROWIN,COLIN);[ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN);
   [ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,ROWIN,COLIN]=reduceAll(ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,ROWIN,COLIN);[ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN);

    S1 = updateTen_CH(MASTER,RDC1,RDC2, VECTORS, VECTORS_CH);
    RP1 = NVR_RDC2PROB(ASSIGNTABLE,RDC1,VECTORS,S1, ROWIN, COLIN);
    RP2 = NVR_RDC2PROB_CH(ASSIGNTABLE,RDC2,VECTORS_CH,S1, ROWIN, COLIN);   
   
   
end





