function [assignmentAccuracy,NVR_SCORE,oldNVR_SCORE]= NVR(peaks,rdcs,HDEXCHANGE, peakIDs, NOES, VECTORS,TYPES, ...
	     RESNUMS,SSTRUCT, HBOND, ALLDISTS, SHIFTS_Filename, SHIFTX_Filename)

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

addpath('~njp/code/JBN-Submission-Snapshot-06-15-07/NVR/General');


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
ROWIN       = 1:size(ASSIGNTABLE,1);
COLIN       = 1:size(ASSIGNTABLE,2);
%This is the master assignment table
MASTER      = ASSIGNTABLE*0;

%prune via amide exchange
HDE = NVR_HD2PROB(ASSIGNTABLE,HDEXCHANGE,HBOND);
ASSIGNTABLE = and(ASSIGNTABLE,HDE);
nlast = sum(sum(ASSIGNTABLE));
for(i=1:100)
   [ASSIGNTABLE] = lockdown(MASTER,ASSIGNTABLE,NOES,ALLDISTS, NTH,ROWIN,COLIN);
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
      [MASTER,ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,ROWIN,COLIN]=combine(SXCP.*CP,SSCP,MASTER,ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,ROWIN,COLIN,NOES,ALLDISTS,NTH,S1,S2,RDC1,RDC2,VECTORS,i);
      [MASTER,ASSIGNTABLE]=updateMASTER(MASTER,ASSIGNTABLE,ROWIN,COLIN);[ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN);
      [ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,ROWIN,COLIN]=reduceAll(ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,ROWIN,COLIN);[ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN);
      [ASSIGNTABLE] = lockdown(MASTER,ASSIGNTABLE,NOES,ALLDISTS, NTH,ROWIN,COLIN);
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
   [ASSIGNTABLE] = lockdown(MASTER,ASSIGNTABLE,NOES,ALLDISTS, NTH,ROWIN,COLIN);
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
   [ASSIGNTABLE] = lockdown(MASTER,ASSIGNTABLE,NOES,ALLDISTS, NTH,ROWIN,COLIN);
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
      [ASSIGNTABLE] = lockdown(MASTER,ASSIGNTABLE,NOES,ALLDISTS, NTH,ROWIN,COLIN);
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
      [ASSIGNTABLE] = lockdown(MASTER,ASSIGNTABLE,NOES,ALLDISTS, NTH,ROWIN,COLIN);
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
      [ASSIGNTABLE] = lockdown(MASTER,ASSIGNTABLE,NOES,ALLDISTS, NTH,ROWIN,COLIN);
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
      [ASSIGNTABLE] = lockdown(MASTER,ASSIGNTABLE,NOES,ALLDISTS, NTH,ROWIN,COLIN);
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


fprintf('computing new NVR_SCORE...\n');

ASSIGNTABLE = MASTER*0+1;
ROWIN=1:size(ASSIGNTABLE,1);
COLIN=1:size(ASSIGNTABLE,2);

NTH=5.5;minval=10e-40;
RP1_A = NVR_RDC2PROB(ASSIGNTABLE,RDC1,VECTORS,S1)+minval;RP1_A =ren(RP1_A);
RP2_A = NVR_RDC2PROB(ASSIGNTABLE,RDC2,VECTORS,S2)+minval;RP2_A =ren(RP2_A);
HDE_A = NVR_HD2PROB(ASSIGNTABLE ,HDEXCHANGE,HBOND)+minval;HDE_A=ren(HDE_A);
CP_A = NVR_CS2PROB(ASSIGNTABLE ,HSHIFTS,NSHIFTS,TYPES,SSTRUCT,NOES,ALLDISTS,NTH,ROWIN,COLIN)+minval;CP_A =ren(CP_A);
SXCP_A = NVR_SHIFTX2PROB(ASSIGNTABLE,HSHIFTS,NSHIFTS,TYPES,SSTRUCT, ...
			NOES,ALLDISTS,NTH,ROWIN,COLIN, SHIFTX_Filename)+minval;SXCP_A =ren(SXCP_A);
SSCP_A = NVR_SHIFTS2PROB(ASSIGNTABLE,HSHIFTS,NSHIFTS,TYPES,SSTRUCT, ...
			NOES,ALLDISTS,NTH,ROWIN,COLIN, SHIFTS_Filename, ...
			SHIFTX_Filename)+minval;SSCP_A = ...
	 ren(SSCP_A);
NVR_SCORE = computeNVR_SCORE(CP_A,SXCP_A,SSCP_A,RP1_A,RP2_A,MASTER)
JOINT =ren(RP1_A.*RP2_A.*HDE_A.*CP_A.*SXCP_A.*SSCP_A).*MASTER;
JOINT = [nonzeros(JOINT)'];
oldNVR_SCORE=mean(log(JOINT))

assignmentAccuracy = computeAssignmentAccuracy(peakIDs, RESNUMS, MASTER);



return






function [ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,ROWIN,COLIN]=reduceAll(ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,ROWIN,COLIN);
[ASSIGNTABLE,CP,x,y]=reduce(ASSIGNTABLE,CP,ROWIN,COLIN);
[ASSIGNTABLE,SXCP,x,y]=reduce(ASSIGNTABLE,SXCP,ROWIN,COLIN);
[ASSIGNTABLE,SSCP,x,y]=reduce(ASSIGNTABLE,SSCP,ROWIN,COLIN);
[ASSIGNTABLE,ASSIGNTABLE,ROWIN,COLIN]=reduce(ASSIGNTABLE,ASSIGNTABLE,ROWIN,COLIN);
ASSIGNTABLE=and(ASSIGNTABLE,ASSIGNTABLE);

function  [MASTER,ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,ROWIN,COLIN]=combine(A,B,MASTER,ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,ROWIN,COLIN,NOES,ALLDISTS,NTH,S1,S2,RDC1,RDC2,VECTORS,THR);
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
   S1 = updateTen(MASTER,S1,RDC1,VECTORS);
   S2 = updateTen(MASTER,S2,RDC2,VECTORS);
   RP1 = NVR_RDC2PROB(ASSIGNTABLE,RDC1,VECTORS,S1);
   RP2 = NVR_RDC2PROB(ASSIGNTABLE,RDC2,VECTORS,S2);
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








function V = computeNVR_SCORE(CP,SXCP,SSCP,RP1,RP2,MASTER)
%singles
[A]=getNVR_SCORE(CP,CP,MASTER);
[B]=getNVR_SCORE(SXCP,SXCP, MASTER);
[C]=getNVR_SCORE(SSCP,SSCP, MASTER);
[D]=getNVR_SCORE(RP1,RP1, MASTER);
[E]=getNVR_SCORE(RP2,RP2, MASTER);
%doubles
[F]=getNVR_SCORE(SXCP.*CP,SXCP.*CP, MASTER);
[G]=getNVR_SCORE(SXCP.*RP1,SXCP.*RP1, MASTER);
[H]=getNVR_SCORE(SXCP.*RP2,SXCP.*RP2, MASTER);
[I]=getNVR_SCORE(SXCP.*SSCP,SXCP.*SSCP, MASTER);
[J]=getNVR_SCORE(SSCP.*CP,SSCP.*CP, MASTER);
[K]=getNVR_SCORE(SSCP.*RP1,SSCP.*RP1, MASTER);
[L]=getNVR_SCORE(SSCP.*RP2,SSCP.*RP2, MASTER);
[M]=getNVR_SCORE(CP.*RP1,CP.*RP1, MASTER);
[N]=getNVR_SCORE(CP.*RP2,CP.*RP2, MASTER);
[O]=getNVR_SCORE(RP1.*RP2,RP1.*RP2, MASTER);
%triples
[P]=getNVR_SCORE(CP.*SXCP.*SSCP,CP.*SXCP.*SSCP, MASTER);
[Q]=getNVR_SCORE(CP.*SXCP.*RP1,CP.*SXCP.*RP1, MASTER);
[R]=getNVR_SCORE(CP.*SXCP.*RP2,CP.*SXCP.*RP2, MASTER);
[S]=getNVR_SCORE(CP.*SSCP.*RP1,CP.*SSCP.*RP1, MASTER);
[T]=getNVR_SCORE(CP.*SSCP.*RP2,CP.*SSCP.*RP2, MASTER);
[U]=getNVR_SCORE(SXCP.*SSCP.*RP1,SXCP.*SSCP.*RP1, MASTER);
[V]=getNVR_SCORE(SXCP.*SSCP.*RP2,SXCP.*SSCP.*RP2, MASTER);
[W]=getNVR_SCORE(SXCP.*RP1.*RP2,SXCP.*RP1.*RP2, MASTER);
[X]=getNVR_SCORE(SSCP.*RP1.*RP2,SSCP.*RP2.*RP1, MASTER);
[Y]=getNVR_SCORE(CP.*RP1.*RP2,CP.*RP2.*RP1, MASTER);
%quads
[Z]=getNVR_SCORE(CP.*RP1.*RP2.*SXCP,CP.*RP2.*RP1.*SXCP, MASTER);
[AA]=getNVR_SCORE(CP.*RP1.*RP2.*SSCP,CP.*RP2.*RP1.*SSCP, MASTER);
[BB]=getNVR_SCORE(SXCP.*RP1.*RP2.*SSCP,SXCP.*RP2.*RP1.*SSCP, MASTER);
[CC]=getNVR_SCORE(SXCP.*CP.*RP2.*SSCP,SXCP.*RP2.*SSCP.*CP, MASTER);
[DD]=getNVR_SCORE(SXCP.*CP.*RP1.*SSCP,SXCP.*RP1.*SSCP.*CP, MASTER);
[EE]=getNVR_SCORE(RP1.*SXCP.*CP.*RP2.*SSCP,RP1.*SXCP.*RP2.*SSCP.*CP, MASTER);
V = A+B+C+D+E+F+G+H+I+J+K+L+M+N+O+P+Q+R+S+T+U+V+W+X+Y+Z+AA+BB+CC+DD+EE;


function score=getNVR_SCORE(bpg,dummy, MASTER);
bpg = ren(bpg) .* MASTER;
bpg = [nonzeros(bpg)'];
score = mean(log(bpg));%sum rather than the mean since the latter
                      %does not consider some assignments
                      %with 0 edge weight in the current bpg.
