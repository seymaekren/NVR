function NVR(peaks,rdcs,HDEXCHANGE, peakIDs, NOES, VECTORS,TYPES,RESNUMS,SSTRUCT, HBOND, ALLDISTS)

%NVR: This program computes resonance assignments from a set of experimental data and a model 
%      Input:  peaks = a Nx2 matrix containing N hsqc peaks. column 1 is the H shift, column 2 is the N shift
%					rdcs = a Nx2 matrix containing N rdcs in 2 media medium, in an arbitrary order.
%					HDEXCHANGE= a Nx1 matrix containing N boolean values indicating whether that peak is a slow-exchanging peak. 1 = slow, 0 = fast
%					peakIDs= a Nx1 matrix containing N ids for the peak. eg peakIDs(1) = the id for peaks(1)
%								NOTE: the rows of peaks, rdcs, HDEXCHANGE, and NOES must have the same order as the
%								rows of peakIDs. 
%					NOES = a NxN matrix. If element NOES(i,j) = 1, then peakIDs(i) is involved in a Dnn with peakIDs(j)
%					VECTORS = A Mx3 matrix containing the normalized backbone amide bond vectors
%					TYPES = A Mx1 matrix containing the 3 letter amino acid code for the backbone amide bond vectors
%					RESNUMS = A Mx1 matrix containing the resiude. Note the rows in the parameters named 
%								 VECTORS, TYPES, SSTRUCT, HBOND, and ALLDISTS should all have the same order
% 								 as listed in the rows of RESNUMS
%					SSTRUCT = A Mx1 matrix containing secondary structure type each residue in the model. C=coil, B=beta, H=helix
%					HBOND = A Mx1 matrix containing the distances between all pairs of backbone amide protons, Y=yes, its involved in an H bond (or is not solvent accesible), N=means that it is essntially labile. 
%					ALLDISTS = A MxM matrix containing the distances between all pairs of backbone amide protons
%		 Output:   


%////////////////////////////////////////////////////////////////////////////////////////////
%//  NVR.m
%//
%//  Version:		0.1
%//
%//  Description:	 This program computes resonance assignments from a set of experimental data and a model
%//
%// authors:
%//    initials    name            organization 					email
%//   ---------   --------------  ------------------------    ------------------------------
%//     CJL         Chris Langmead  Dartmouth College         langmead@dartmouth.edu
%//
%//
%// history:
%//     when        who     what
%//     --------    ----    ----------------------------------------------------------
%//     12/02/03    CJL 	 initial version for publication [Langmead et al, J Biomol NMR 2004]
%//
%////////////////////////////////////////////////////////////////////////////////////////////

%    NVR
%    This library is free software; you can redistribute it and/or
%    modify it under the terms of the GNU Lesser General Public
%    License as published by the Free Software Foundation; either
%    version 2.1 of the License, or (at your option) any later version.

%    This library is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%    Lesser General Public License for more details.

%    You should have received a copy of the GNU Lesser General Public
%    License along with this library; if not, write to the Free Software
%    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

% 		Contact Info:
%							Bruce Randall Donald
%							HB 6211
%							Dartmouth College
%							Hanover, NH 03755
%							brd@cs.dartmouth.edu

% 		If you use publish any results derived from the use of this program please cite:
%		"An Expectation/Maximization Nuclear Vector Replacement Algorithm for Automated NMR Resonance Assignments," 
%		C. J. Langmead and B. R. Donald, 
%		Journal of Biomolecular NMR, 2004 (in press)


%  Copyright (C) 2003  Christopher James Langmead and Bruce R. Donald
%
%  <signature of Bruce Donald>, 2 December 2003
%  Bruce Donald, Professor of Computer Science

fprintf('NOTE: This script is neither fully tested nor intended for general use as a piece of software. ');
fprintf('It is only guaranteed to work on the example data provided. You are encouraged to modify and ');
fprintf('improve it for your own needs. The original authors are presently implementing more robust and ');
fprintf('general purpose software based on the algorithms presented in this script.\n');



HSHIFTS = peaks(:,1);
NSHIFTS = peaks(:,2);
RDC1 = rdcs(:,1);
RDC2 = rdcs(:,2);

%compute a tollerance for NOE distances
NTH=4.8;
mu=mean(mean(ALLDISTS));
if(mu-12.9>0)
   NTH=NTH+(mu-12.9); NTH=min(NTH,8);
end

ASSIGNTABLE = ones(length(HSHIFTS),size(VECTORS,1))/size(VECTORS,1);
%these keep track of which peaks and residues are represented in the current
%matricies
ROWIN=1:size(ASSIGNTABLE,1);
COLIN=1:size(ASSIGNTABLE,2);
%This is the master assignment table
MASTER=ASSIGNTABLE*0;

%prune via amide exchange
HDE = NVR_HD2PROB(ASSIGNTABLE,HDEXCHANGE,HBOND);
ASSIGNTABLE = and(ASSIGNTABLE,HDE);
nlast = sum(sum(ASSIGNTABLE));
for(i=1:100)
   [ASSIGNTABLE] = lockdown(ASSIGNTABLE,NOES,ALLDISTS, NTH,ROWIN,COLIN);
   if(sum(sum(ASSIGNTABLE )) == nlast)
      break;
   end
   nlast = sum(sum(ASSIGNTABLE ));
end

%prune via BMRB stats
[CP] = NVR_CS2PROB(ASSIGNTABLE,HSHIFTS,NSHIFTS,TYPES,SSTRUCT,NOES,ALLDISTS,NTH,ROWIN,COLIN);SXCP=CP;SSCP = CP;RP1 = CP;RP2 = CP;
ASSIGNTABLE = and(ASSIGNTABLE,CP);
for(lv = [round(size(CP,2)*.6)])
   V = vote(CP,CP,CP,CP,CP,ASSIGNTABLE,NOES,ALLDISTS,NTH,ROWIN,COLIN,2,lv);
   V2 = vote(CP',CP',CP',CP',CP',ASSIGNTABLE',NOES,ALLDISTS,NTH,ROWIN,COLIN,2,lv);
   V = or(V,V2');
   ASSIGNTABLE = ASSIGNTABLE.*V(1:size(ASSIGNTABLE,1),1:size(ASSIGNTABLE,2));
   
   nlast = sum(sum(ASSIGNTABLE));
   for(i=1:100)
      [ASSIGNTABLE] = lockdown(ASSIGNTABLE,NOES,ALLDISTS, NTH,ROWIN,COLIN);
      if(sum(sum(ASSIGNTABLE)) == nlast)
         break;
      end
      nlast = sum(sum(ASSIGNTABLE));
   end
end

%prune via the program shiftx
[SXCP] = NVR_SHIFTX2PROB(ASSIGNTABLE,HSHIFTS,NSHIFTS,TYPES,SSTRUCT,NOES,ALLDISTS,NTH,ROWIN,COLIN);
ASSIGNTABLE = and(ASSIGNTABLE,SXCP);
for(lv = [round(size(CP,2)*.5)])
   %reduce to a fixed percentage
   V = vote(SXCP,SXCP,SXCP,SXCP,SXCP,ASSIGNTABLE,NOES,ALLDISTS,NTH,ROWIN,COLIN,2,lv);
   V2 = vote(SXCP',SXCP',SXCP',SXCP',SXCP',ASSIGNTABLE',NOES,ALLDISTS,NTH,ROWIN,COLIN,2,lv);
   V = or(V,V2');
   ASSIGNTABLE = ASSIGNTABLE.*V(1:size(ASSIGNTABLE,1),1:size(ASSIGNTABLE,2));
   nlast = sum(sum(ASSIGNTABLE));
   for(i=1:100)
      [ASSIGNTABLE] = lockdown(ASSIGNTABLE,NOES,ALLDISTS, NTH,ROWIN,COLIN);
      if(sum(sum(ASSIGNTABLE)) == nlast)
         break;
      end
      
      nlast = sum(sum(ASSIGNTABLE));
   end
end



%prune via the program shifts
[SSCP] = NVR_SHIFTS2PROB(ASSIGNTABLE,HSHIFTS,NSHIFTS,TYPES,SSTRUCT,NOES,ALLDISTS,NTH,ROWIN,COLIN);
ASSIGNTABLE = and(ASSIGNTABLE,SSCP);


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
      [ASSIGNTABLE] = lockdown(ASSIGNTABLE,NOES,ALLDISTS, NTH,ROWIN,COLIN);
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
[ASSIGNTABLE] = lockdown(ASSIGNTABLE,NOES,ALLDISTS, NTH,ROWIN,COLIN);
[MASTER,ASSIGNTABLE]=updateMASTER(MASTER,ASSIGNTABLE,ROWIN,COLIN);[ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN);
[ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,ROWIN,COLIN]=reduceAll(ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,ROWIN,COLIN);[ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN);

if(sum(sum(MASTER))>=5)
   S1 = ones(3);S2=S1;
   S1 = updateTen(MASTER,S1,RDC1,VECTORS);
   S2 = updateTen(MASTER,S2,RDC2,VECTORS);
   RP1 = NVR_RDC2PROB(ASSIGNTABLE,RDC1,VECTORS,S1);
   RP2 = NVR_RDC2PROB(ASSIGNTABLE,RDC2,VECTORS,S2);
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
      [ASSIGNTABLE] = lockdown(ASSIGNTABLE,NOES,ALLDISTS, NTH,ROWIN,COLIN);
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
      [MASTER,ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,ROWIN,COLIN]=combine(SXCP.*CP,SSCP,MASTER,ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,ROWIN,COLIN,NOES,ALLDISTS,NTH,S1,S2,RDC1,RDC2,VECTORS,i);
      [MASTER,ASSIGNTABLE]=updateMASTER(MASTER,ASSIGNTABLE,ROWIN,COLIN);[ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN);
      [ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,ROWIN,COLIN]=reduceAll(ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,ROWIN,COLIN);[ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN);
      [ASSIGNTABLE] = lockdown(ASSIGNTABLE,NOES,ALLDISTS, NTH,ROWIN,COLIN);
      [MASTER,ASSIGNTABLE]=updateMASTER(MASTER,ASSIGNTABLE,ROWIN,COLIN);[ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN);
      [ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,ROWIN,COLIN]=reduceAll(ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,ROWIN,COLIN);[ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN);
      
      S1 = updateTen(MASTER,S1,RDC1,VECTORS);
      S2 = updateTen(MASTER,S2,RDC2,VECTORS);
      RP1 = NVR_RDC2PROB(ASSIGNTABLE,RDC1,VECTORS,S1);
      RP2 = NVR_RDC2PROB(ASSIGNTABLE,RDC2,VECTORS,S2);
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
   [MASTER,ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,ROWIN,COLIN]=combine(SXCP.*CP,SSCP,MASTER,ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,ROWIN,COLIN,NOES,ALLDISTS,NTH,S1,S2,RDC1,RDC2,VECTORS,.9);
   [MASTER,ASSIGNTABLE]=updateMASTER(MASTER,ASSIGNTABLE,ROWIN,COLIN);[ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN);
   [ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,ROWIN,COLIN]=reduceAll(ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,ROWIN,COLIN);[ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN);
   [ASSIGNTABLE] = lockdown(ASSIGNTABLE,NOES,ALLDISTS, NTH,ROWIN,COLIN);
   [MASTER,ASSIGNTABLE]=updateMASTER(MASTER,ASSIGNTABLE,ROWIN,COLIN);[ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN);
   [ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,ROWIN,COLIN]=reduceAll(ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,ROWIN,COLIN);[ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN);
   
   S1 = updateTen(MASTER,S1,RDC1,VECTORS);
   S2 = updateTen(MASTER,S2,RDC2,VECTORS);
   RP1 = NVR_RDC2PROB(ASSIGNTABLE,RDC1,VECTORS,S1);
   RP2 = NVR_RDC2PROB(ASSIGNTABLE,RDC2,VECTORS,S2);
   if(size(MASTER,1)-sum(sum(MASTER))<=1)
      break;
   end
   if(sum(sum(ASSIGNTABLE))==ilast)
      break;
   end
   ilast = sum(sum(ASSIGNTABLE));
end

S1 = updateTen(MASTER,S1,RDC1,VECTORS);
S2 = updateTen(MASTER,S2,RDC2,VECTORS);
RP1 = NVR_RDC2PROB(ASSIGNTABLE,RDC1,VECTORS,S1);
RP2 = NVR_RDC2PROB(ASSIGNTABLE,RDC2,VECTORS,S2);
for(q=[.9999,.999,.99,.9])
   [MASTER,ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,ROWIN,COLIN]=combine(RP1,RP2,MASTER,ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,ROWIN,COLIN,NOES,ALLDISTS,NTH,S1,S2,RDC1,RDC2,VECTORS,q);[ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN);
   [MASTER,ASSIGNTABLE]=updateMASTER(MASTER,ASSIGNTABLE,ROWIN,COLIN);[ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN);
   [ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,ROWIN,COLIN]=reduceAll(ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,ROWIN,COLIN);[ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN);
   [ASSIGNTABLE] = lockdown(ASSIGNTABLE,NOES,ALLDISTS, NTH,ROWIN,COLIN);
   [MASTER,ASSIGNTABLE]=updateMASTER(MASTER,ASSIGNTABLE,ROWIN,COLIN);[ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN);
   [ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,ROWIN,COLIN]=reduceAll(ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,ROWIN,COLIN);[ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN);
   
   S1 = updateTen(MASTER,S1,RDC1,VECTORS);
   S2 = updateTen(MASTER,S2,RDC2,VECTORS);
   RP1 = NVR_RDC2PROB(ASSIGNTABLE,RDC1,VECTORS,S1);
   RP2 = NVR_RDC2PROB(ASSIGNTABLE,RDC2,VECTORS,S2);
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
      [ASSIGNTABLE] = lockdown(ASSIGNTABLE,NOES,ALLDISTS, NTH,ROWIN,COLIN);
      [MASTER,ASSIGNTABLE]=updateMASTER(MASTER,ASSIGNTABLE,ROWIN,COLIN);[ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN);
      [ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,ROWIN,COLIN]=reduceAll(ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,ROWIN,COLIN);[ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN);
      S1 = updateTen(MASTER,S1,RDC1,VECTORS);
      S2 = updateTen(MASTER,S2,RDC2,VECTORS);
      RP1 = NVR_RDC2PROB(ASSIGNTABLE,RDC1,VECTORS,S1);
      RP2 = NVR_RDC2PROB(ASSIGNTABLE,RDC2,VECTORS,S2);
      if(sum(sum(ASSIGNTABLE)) == nlast)
         break;
      end
      nlast = sum(sum(ASSIGNTABLE));
   end
end

%phase 2
lv = 5;
for(q=1:size(MASTER,1))
   if(sum(sum(MASTER))==size(MASTER,1))
      break;
   end
   
   %do voting
   V = vote(CP,SXCP,SSCP,RP1,RP2,ASSIGNTABLE,NOES,ALLDISTS,NTH,ROWIN,COLIN,1,0);oV=V;
   
   %eliminate any that we know can't be possible
   V(1:size(ASSIGNTABLE,1),1:size(ASSIGNTABLE,2))=V(1:size(ASSIGNTABLE,1),1:size(ASSIGNTABLE,2)).*ASSIGNTABLE;
   
   %compute the assignment, given the distribution in the vote table
   [A]=getASS(V,V, ASSIGNTABLE,1,1);A(1:size(ASSIGNTABLE,1),1:size(ASSIGNTABLE,2))=A(1:size(ASSIGNTABLE,1),1:size(ASSIGNTABLE,2)).*ASSIGNTABLE;
   %only consider the ones with the max num votes
   V = thresh(oV,max(max(max(V)),1));
   
   %threshold anything that is really unlikely
   nV=RP1.*RP2.*CP.*SXCP.*SSCP;
   nV=thresh(nV,min(nonzeros(nV))*1);
   dontrun=0;
   
   if(max(max(V))<30 & size(ASSIGNTABLE,1)>3)
      nV=thresh(nV,mean(nonzeros(nV)));
   end
   
   %take intersection of A and V and nV
   X = and(A(1:size(ASSIGNTABLE,1),1:size(ASSIGNTABLE,2)),and(V(1:size(ASSIGNTABLE,1),1:size(ASSIGNTABLE,2)),nV));
   if(sum(sum(X))==0)
      ct = 1;
      while(sum(sum(X)) < 1 & ct<31)
         V = thresh(oV,max(max(max(V))-ct,1));ct=ct+1;
         X = and(A(1:size(ASSIGNTABLE,1),1:size(ASSIGNTABLE,2)),and(V(1:size(ASSIGNTABLE,1),1:size(ASSIGNTABLE,2)),nV));
      end
   end
   
   %make the assignments
   for(i=1:size(X,1))
      x = find(X(i,:));
      if(length(x)==1)
         ASSIGNTABLE(i,:)=0;
         ASSIGNTABLE(i,x)=1;
      end
   end
   
   
   nlast = sum(sum(ASSIGNTABLE));
   for(i=1:100)
      [MASTER,ASSIGNTABLE]=updateMASTER(MASTER,ASSIGNTABLE,ROWIN,COLIN);[ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN);
      [ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,ROWIN,COLIN]=reduceAll(ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,ROWIN,COLIN);[ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN);
      [ASSIGNTABLE] = lockdown(ASSIGNTABLE,NOES,ALLDISTS, NTH,ROWIN,COLIN);
      [MASTER,ASSIGNTABLE]=updateMASTER(MASTER,ASSIGNTABLE,ROWIN,COLIN);[ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN);
      [ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,ROWIN,COLIN]=reduceAll(ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,ROWIN,COLIN);[ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN);
      S1 = updateTen(MASTER,S1,RDC1,VECTORS);
      S2 = updateTen(MASTER,S2,RDC2,VECTORS);
      RP1 = NVR_RDC2PROB(ASSIGNTABLE,RDC1,VECTORS,S1);
      RP2 = NVR_RDC2PROB(ASSIGNTABLE,RDC2,VECTORS,S2);
      if(sum(sum(ASSIGNTABLE)) == nlast)
         break;
      end
      nlast = sum(sum(ASSIGNTABLE));
   end
   if(size(MASTER,1)-sum(sum(MASTER))<=1)
      break;
   end
   if(sum(sum(ASSIGNTABLE))==size(ASSIGNTABLE,1))
      break;
   end
   
   for(qq=[1,.8,1,.8,1,.8,1,.8])
      mx = max(max(RP1.*RP2));
      [MASTER,ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,ROWIN,COLIN]=combine(RP1,RP2,MASTER,ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,ROWIN,COLIN,NOES,ALLDISTS,NTH,S1,S2,RDC1,RDC2,VECTORS,qq*mx);[ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN);
      
      [MASTER,ASSIGNTABLE]=updateMASTER(MASTER,ASSIGNTABLE,ROWIN,COLIN);[ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN);
      [ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,ROWIN,COLIN]=reduceAll(ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,ROWIN,COLIN);[ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN);
      [ASSIGNTABLE] = lockdown(ASSIGNTABLE,NOES,ALLDISTS, NTH,ROWIN,COLIN);
      [MASTER,ASSIGNTABLE]=updateMASTER(MASTER,ASSIGNTABLE,ROWIN,COLIN);[ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN);
      [ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,ROWIN,COLIN]=reduceAll(ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,ROWIN,COLIN);[ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN);
      
      S1 = updateTen(MASTER,S1,RDC1,VECTORS);
      S2 = updateTen(MASTER,S2,RDC2,VECTORS);
      %[Rotations(S1,R1) Rotations(S2,R2)]
      RP1 = NVR_RDC2PROB(ASSIGNTABLE,RDC1,VECTORS,S1);
      RP2 = NVR_RDC2PROB(ASSIGNTABLE,RDC2,VECTORS,S2);
      if(size(MASTER,1)-sum(sum(MASTER))<=1)
         break;
      end
   end
   
   if(size(MASTER,1)-sum(sum(MASTER))<=1)
      break;
   end
   
   %reduce again
   V = vote(CP,SXCP,SSCP,RP1,RP2,ASSIGNTABLE,NOES,ALLDISTS,NTH,ROWIN,COLIN,2,lv);
   V2 = vote(CP',SXCP',SSCP',RP1',RP2',ASSIGNTABLE',NOES,ALLDISTS,NTH,ROWIN,COLIN,2,lv);
   lv = min(lv-4,1);
   V = or(V,V2');
   ASSIGNTABLE = ASSIGNTABLE.*V(1:size(ASSIGNTABLE,1),1:size(ASSIGNTABLE,2));
   nlast = sum(sum(ASSIGNTABLE));
   for(i=1:100)
      [MASTER,ASSIGNTABLE]=updateMASTER(MASTER,ASSIGNTABLE,ROWIN,COLIN);[ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN);
      [ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,ROWIN,COLIN]=reduceAll(ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,ROWIN,COLIN);[ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN);
      [ASSIGNTABLE] = lockdown(ASSIGNTABLE,NOES,ALLDISTS, NTH,ROWIN,COLIN);
      [MASTER,ASSIGNTABLE]=updateMASTER(MASTER,ASSIGNTABLE,ROWIN,COLIN);[ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN);
      [ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,ROWIN,COLIN]=reduceAll(ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,ROWIN,COLIN);[ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN);
      S1 = updateTen(MASTER,S1,RDC1,VECTORS);
      S2 = updateTen(MASTER,S2,RDC2,VECTORS);
      RP1 = NVR_RDC2PROB(ASSIGNTABLE,RDC1,VECTORS,S1);
      RP2 = NVR_RDC2PROB(ASSIGNTABLE,RDC2,VECTORS,S2);
      if(sum(sum(ASSIGNTABLE)) == nlast)
         break;
      end
      nlast = sum(sum(ASSIGNTABLE));
   end
end

if(size(MASTER,1)-sum(sum(MASTER))==1)
   %make last assignment
   for(i=1:size(MASTER,1))
      x = find(MASTER(i,:));
      if(length(x)==0)
         mx=i;
         break;
      end
   end
   for(i=1:size(MASTER,2))
      x = find(MASTER(:,i));
      if(length(x)==0)
         my=i;
         break;
      end
   end
   MASTER(mx,my)=1;
end

%now, compute correctness

a = load('answerkey.m');
assignments=0;
correct=0;
incorr=0;

for(i=1:size(MASTER,1))
   pk=peakIDs(i);%get the peak id
   x = find(MASTER(i,:));
    
   rn=RESNUMS(x);%get the residue it was assigned to 
   
   assignments(i,1)=pk;
   assignments(i,2)=rn;
  
   foo = find(a(:,1)==pk);
   if(rn==a(foo,2))
      correct=correct+1;
      incorr(i)=0;
    else
      incorr(i)=1;
   end

end



fprintf('\n');
fprintf('\n');
fprintf('Assignment Accuracy = %f percent \n',correct/size(MASTER,1)*100);
fprintf('\n');
fprintf('\n');
fprintf('Assignments\n');
fprintf('Peak ID -> Residue Number\n');
for(i=1:size(assignments,1))
   if(incorr(i)==0)
      fprintf('%d	%d\n',assignments(i,1),assignments(i,2));
   else
      fprintf('*%d	%d  (correct=%d->%d)\n',assignments(i,1),assignments(i,2),a(i,1),a(i,2));
   end
end
fprintf('\n');
fprintf('\n');

if(size(MASTER,1)~=size(MASTER,2))
   %find which residues aren't assigned
   ct=1;
   missing=[];
   
   for(i=1:length(RESNUMS))
      rn = RESNUMS(i);
      if(length(find(assignments(:,2)==rn))==0)
         missing(ct)=rn;
         ct=ct+1;
      end
   end
   
   fprintf('The following residues are missing peaks\n');
   for(i=1:length(missing))
     fprintf('%d\n',missing(i));
   end
end





return


function [ASSIGNTABLE]=getUACMB(A,B, ASSIGNTABLE,NOES,ALLDISTS,NTH,THR,rep,ROWIN,COLIN);
tot=0;
for(k=1:rep)
   %find unambiguous
   for(i=1:size(A,1))
      A(i,:)=A(i,:)/sum(A(i,:));
      B(i,:)=B(i,:)/sum(B(i,:));
   end
   UA = NVR_FINDUACMB(A,B,ASSIGNTABLE,THR);
   if(sum(sum(UA))==tot)
      break;
   else
      tot = sum(sum(UA));
   end
   %make unambiguous assignments
   for(i=1:size(UA,1))
      x = find(UA(i,:));
      if(length(x)>0)
         ASSIGNTABLE(i,:)=0;
         ASSIGNTABLE(:,x)=0;
         ASSIGNTABLE(i,x)=1;
      end
   end
   %apply noes
   [ASSIGNTABLE] = lockdown(ASSIGNTABLE,NOES,ALLDISTS,NTH,ROWIN,COLIN);
end


function [ASSIGNTABLE] = lockdown(ASSIGNTABLE,NOES,ALLDISTS, NTH,ROWIN,COLIN)
%prune via NOEs
NP = NVR_NOE2PROB(ASSIGNTABLE,NOES,ALLDISTS,NTH,ROWIN,COLIN);
ASSIGNTABLE=and(NP,NP);


%##################################################################################
%##  updateTen
function SAUPE = updateTen(TABLE,ET,RDC,VECS)

SAUPE  = ET;

M = zeros(1,5);
vec = 0;
ct=1;
for(i=1:size(TABLE,1))
   x = find(TABLE(i,:));
   if(length(x)==1 & RDC(i)>-999)
      M(ct,1) = VECS(x,1).^2-VECS(x,3).^2;
      M(ct,2) = 2*(VECS(x,1)*VECS(x,2));
      M(ct,3) = 2*(VECS(x,1)*VECS(x,3));
      M(ct,4) = VECS(x,2).^2-VECS(x,3).^2;
      M(ct,5) = 2*(VECS(x,2)*VECS(x,3));
      vec(ct) =  RDC(i); 
      ct=ct+1;
   end
end
if(size(M,1)<5)
   return
end
%invert the matrix
if(size(M,1)==5 & size(M,2)==5)
   M = inv(M);   
else
   M = pinv(M);
end
%do the least squares fitting
s = M*vec';
SAUPE = zeros(3,3);
SAUPE(1,1) = s(1);
SAUPE(2,2) = s(4);
SAUPE(3,3) = -1*(s(1)+s(4));
SAUPE(1,2) = s(2);
SAUPE(2,1) = s(2);
SAUPE(1,3) = s(3);
SAUPE(3,1) = s(3);
SAUPE(2,3) = s(5);
SAUPE(3,2) = s(5);





function [MASTER,A]=updateMASTER(MASTER,A,ROWIN,COLIN);
%first, propagate contraints
recur=1;
rm=ones(1,size(A,1));
cm=ones(1,size(A,2));
ain=A;
while(recur==1)
   recur=0;
   for(i=1:size(A,1))
      x=find(A(i,:));
      if(length(x)==1)
         A(:,x)=0;
         A(i,x)=1;
         if(rm(i)==1)
            rm(i)=0;
            recur=1;
         end
      end
   end
   for(i=1:size(A,2))
      x=find(A(:,i));
      if(length(x)==1 & size(A,1)==size(A,2))
         A(x,:)=0;
         A(x,i)=1;
         if(cm(i)==1)
            cm(i)=0;
            recur=1;
         end
      end
   end
end
for(i=1:size(A,1))
   if(sum(A(i,:))==0)
      A(i,:)=ain(i,:);
      x=find(A(i,:));
      A(x,:)=ain(x,:);
   end
end

%next, make assignments
for(i=1:size(A,1))
   x = find(A(i,:));
   if(length(x)==1)
      MASTER(ROWIN(i),COLIN(x))=1;
   end
end


function [ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,ROWIN,COLIN]=reduceAll(ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,ROWIN,COLIN);
[ASSIGNTABLE,CP,x,y]=reduce(ASSIGNTABLE,CP,ROWIN,COLIN);
[ASSIGNTABLE,SXCP,x,y]=reduce(ASSIGNTABLE,SXCP,ROWIN,COLIN);
[ASSIGNTABLE,SSCP,x,y]=reduce(ASSIGNTABLE,SSCP,ROWIN,COLIN);
[ASSIGNTABLE,ASSIGNTABLE,ROWIN,COLIN]=reduce(ASSIGNTABLE,ASSIGNTABLE,ROWIN,COLIN);
ASSIGNTABLE=and(ASSIGNTABLE,ASSIGNTABLE);

function [A,B,ROWIN,COLIN]=reduce(A,B,ROWIN,COLIN)
ain=A;
bin=B;
ct=1;row=[1:size(A,1)];col=[1:size(A,2)];
for(i=1:size(A,1))
   x=find(A(i,:));
   if(length(x)==1)
      row(i)=0;
      col(x)=0;   
   end
end
B=B.*A;
for(i=1:size(B,1))
   if(sum(B(i,:))==0)
      B(i,:)=bin(i,:);
   end
   if(sum(B(i,:))==0)
      B(i,:)=1;
   end
   B(i,:)=B(i,:)/sum(B(i,:));
end
if(isempty(row))
   return
end
row=find(row);
col=find(col);
B=B(row,col);
for(i=1:size(B,1))
   if(sum(B(i,:))==0)
      B(i,:)=bin(i,:);
   end
   B(i,:)=B(i,:)/sum(B(i,:));
end
ROWIN=ROWIN(row);
COLIN=COLIN(col);


function  [MASTER,ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,ROWIN,COLIN]=combine(A,B,MASTER,ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,ROWIN,COLIN,NOES,ALLDISTS,NTH,S1,S2,RDC1,RDC2,VECTORS,THR);
THR;
if(THR>0)
   [ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN);
   [ASSIGNTABLE] = getUACMB(A,B,ASSIGNTABLE,NOES,ALLDISTS,NTH,THR,1,ROWIN,COLIN);[ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN);
   [MASTER,ASSIGNTABLE]=updateMASTER(MASTER,ASSIGNTABLE,ROWIN,COLIN);[ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN);
   [ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,ROWIN,COLIN]=reduceAll(ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,ROWIN,COLIN);[ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN);
   [ASSIGNTABLE] = lockdown(ASSIGNTABLE,NOES,ALLDISTS, NTH,ROWIN,COLIN);
   [MASTER,ASSIGNTABLE]=updateMASTER(MASTER,ASSIGNTABLE,ROWIN,COLIN);[ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN);
   [ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,ROWIN,COLIN]=reduceAll(ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,ROWIN,COLIN);[ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN);
   S1 = updateTen(MASTER,S1,RDC1,VECTORS);
   S2 = updateTen(MASTER,S2,RDC2,VECTORS);
   RP1 = NVR_RDC2PROB(ASSIGNTABLE,RDC1,VECTORS,S1);
   RP2 = NVR_RDC2PROB(ASSIGNTABLE,RDC2,VECTORS,S2);
end

function [ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN)
sct=1;
for(i=1:size(MASTER,1))
   x = find(MASTER(i,:));
   if(length(x)==0)
      ROWIN(sct)=i;
      sct=sct+1;
   end
end
sct=1;
for(i=1:size(MASTER,2))
   x = find(MASTER(:,i));
   if(length(x)==0)
      COLIN(sct)=i;
      sct=sct+1;
   end
end


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

function M = getASS(X,Y,ASSIGNTABLE,THR,atype,top)
if(size(ASSIGNTABLE,1)==0)
   M=[];
   return
end
if(size(ASSIGNTABLE,1)-size(ASSIGNTABLE,2)~=0)
   F=ones(max(size(ASSIGNTABLE)));
   F(1:size(ASSIGNTABLE,1),1:size(ASSIGNTABLE,2))=ASSIGNTABLE;
   ASSIGNTABLE=F;   
   for(i=1:size(ASSIGNTABLE,1))
      ASSIGNTABLE(i,:)=ASSIGNTABLE(i,:)/sum(ASSIGNTABLE(i,:));
   end
end
ASSIGNTABLE=and(ASSIGNTABLE,ASSIGNTABLE);
if(size(X,1)-size(X,2)~=0)
   F=ones(max(size(X)));
   F(1:size(X,1),1:size(X,2))=X;
   HM=F;   
   for(i=1:size(HM,1))
      HM(i,:)=HM(i,:)/sum(HM(i,:));
   end
else
   HM=X;
end
M=HM*0;
h=1;
if(atype==1)
   h=hungarian(HM*-1);
elseif(atype==2)
   sz=top;
   for(i=1:size(HM,1))
      [x ind]= sort(HM(i,:));
      in2=find(x);
      ind=ind(in2);x=x(in2);
      if(range(x)>0)
         if(length(x)<sz+1)
            M(i,ind)=1;
         else
            M(i,ind(length(ind)-sz:length(ind)))=1;
         end
      else
         M(i,ind)=1;
      end
      h=[];   
   end
else
   h=simp(HM);
end
for(i=1:length(h))
   M(i,h(i)) = 1;
end

function h = simp(M)
for(i=1:size(M,1))
   x = find(M(i,:));
   if(range(x)>0)
      [m,pos]=max(M(i,:));
      h(i)=pos;
   else
      r =randperm(length(x));
      h(i)=r(1);
   end
end
%throw away any that have more than one peak voting for it
for(i=1:size(M,2))
   x = find(M(:,i));
   if(length(x)>1)
      M(:,i)=0;
   end
end






