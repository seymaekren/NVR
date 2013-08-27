function [HD_SCORE,assignmentAccuracy, assignments, weightedHD_Score] = HD_withTwoVectors(peaks,rdcs,HDEXCHANGE, ...
						  peakIDs, NOES, ...
						  VECTORS,origVECTORS, TYPES, ...
						  RESNUMS,SSTRUCT, HBOND, ALLDISTS, ...
						  IALLDISTS, SHIFTS_Filename, SHIFTX_Filename)

%HD: This program computes how well a given model fits  
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
%					IALLDISTS = A MxM matrix containing an altered version of the distances, that effectively ignores residues
%									involed in random coil
%		 Output:   


%////////////////////////////////////////////////////////////////////////////////////////////
%//  HD.m
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

%    HD
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

fprintf('This is v 1.51\n');

addpath('/home/home4/apaydin/Mist/Matlab/General');
addpath('/home/home4/apaydin/Mist/NVR/Routines');

HSHIFTS = peaks(:,1);
NSHIFTS = peaks(:,2);
RDC1    = rdcs(:,1);
RDC2    = rdcs(:,2);

%compute a tollerance for NOE distances
NTH    = 4.8;
mu     = mean(mean(ALLDISTS));

if(mu-12.9>0)
   NTH = NTH + (mu-12.9); NTH = min(NTH,8);
end

ASSIGNTABLE  = ones(length(HSHIFTS),size(VECTORS,1))/size(VECTORS,1);
OASSIGNTABLE = ASSIGNTABLE;
%these keep track of which peaks and residues are represented in the current
%matricies
ROWIN        = 1:size(ASSIGNTABLE,1);
COLIN        = 1:size(ASSIGNTABLE,2);
%This is the master assignment table
MASTER       = ASSIGNTABLE*0;

fprintf('computing assignments...\n');
HDE          = HD_HD2PROB(ASSIGNTABLE,HDEXCHANGE,HBOND);
ASSIGNTABLE  = and(ASSIGNTABLE,HDE);

TP           = ASSIGNTABLE;
TP           = NVR_TOCSY2PROB(peakIDs,HSHIFTS,NSHIFTS,TYPES,SSTRUCT,NOES,IALLDISTS,NTH,ROWIN,COLIN);
ASSIGNTABLE  = and(ASSIGNTABLE,TP );

[CP]         = HD_CS2PROB(OASSIGNTABLE,HSHIFTS,NSHIFTS,TYPES,SSTRUCT,NOES,IALLDISTS,NTH,ROWIN,COLIN);
ASSIGNTABLE  = and(ASSIGNTABLE,CP);

%prune via the program shiftx
[SXCP]       = HD_SHIFTX2PROB(OASSIGNTABLE,HSHIFTS,NSHIFTS,TYPES,SSTRUCT,NOES,IALLDISTS,NTH,ROWIN,COLIN,SHIFTX_Filename);
ASSIGNTABLE  = and(ASSIGNTABLE,SXCP);   

%prune via the program shifts
[SSCP]       = HD_SHIFTS2PROB(OASSIGNTABLE,HSHIFTS,NSHIFTS,TYPES,SSTRUCT, ...
			NOES,IALLDISTS,NTH,ROWIN,COLIN, SHIFTS_Filename,SHIFTX_Filename);
ASSIGNTABLE  = and(ASSIGNTABLE,SSCP);

SSCP=SSCP.*ASSIGNTABLE;SXCP=SXCP.*ASSIGNTABLE;CP=CP.*ASSIGNTABLE;TP=TP.*ASSIGNTABLE;HDE=HDE.*ASSIGNTABLE;

%initialize the RDC variables
RDC1_VECTORS = VECTORS;
RDC2_VECTORS = origVECTORS;

S1 = updateTen(MASTER,RDC1,RDC1_VECTORS);
S2 = updateTen(MASTER,RDC2,RDC2_VECTORS);
RP1 = NVR_RDC2PROB(ASSIGNTABLE,RDC1,RDC1_VECTORS,S1);
RP2 = NVR_RDC2PROB(ASSIGNTABLE,RDC2,RDC2_VECTORS,S2);

nlast = sum(sum(ASSIGNTABLE));
for(i=1:100)
   [ASSIGNTABLE] = lockdown(MASTER,ASSIGNTABLE,NOES,IALLDISTS, NTH, ROWIN,COLIN);
   [ASSIGNTABLE] = lockdown(MASTER,ASSIGNTABLE,NOES,ALLDISTS, 40, ROWIN,COLIN);
   [MASTER,ASSIGNTABLE]=updateMASTER(MASTER,ASSIGNTABLE,ROWIN,COLIN);[ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN);
   [ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,TP,HDE,ROWIN,COLIN]=reduceAll(ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,TP,HDE,ROWIN,COLIN);[ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN);
   if(sum(sum(ASSIGNTABLE)) == nlast)
      break;
   end
   nlast = sum(sum(ASSIGNTABLE));
end



ilast = sum(sum(ASSIGNTABLE));
for(ww=1:10)
   V = vote5(CP,TP,SSCP,SXCP,HDE,ASSIGNTABLE,NOES,ALLDISTS,NTH,ROWIN,COLIN,1,1);V=thresh(V,1);
   V2 = vote5(1-CP,1-TP,1-SSCP,1-SXCP,1-HDE,ASSIGNTABLE,NOES,ALLDISTS,NTH,ROWIN,COLIN,2,10);
   V3 = vote5(CP,TP,SSCP,SXCP,HDE,ASSIGNTABLE,NOES,ALLDISTS,NTH,ROWIN,COLIN,2,20);
   FOO = or(V2,or(V,V3));ASSIGNTABLE=ASSIGNTABLE.*FOO;
   nlast = sum(sum(ASSIGNTABLE));
   for(i=1:100)
      [ASSIGNTABLE] = lockdown(MASTER,ASSIGNTABLE,NOES,IALLDISTS, NTH, ROWIN,COLIN);
      [ASSIGNTABLE] = lockdown(MASTER,ASSIGNTABLE,NOES,ALLDISTS, 15, ROWIN,COLIN);
      [MASTER,ASSIGNTABLE]=updateMASTER(MASTER,ASSIGNTABLE,ROWIN,COLIN);[ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN);
      [ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,TP,HDE,ROWIN,COLIN]=reduceAll(ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,TP,HDE,ROWIN,COLIN);[ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN);
      if(sum(sum(ASSIGNTABLE)) == nlast)
         break;
      end
      nlast = sum(sum(ASSIGNTABLE));
   end
   if(ilast ==sum(sum(ASSIGNTABLE)))
      break;
   end
   ilast = sum(sum(ASSIGNTABLE));
end

for(lv = [20])
   V = vote5(CP,SXCP,SSCP,TP,HDE,ASSIGNTABLE,NOES,ALLDISTS,NTH,ROWIN,COLIN,2,lv);
   V2 = vote5(1-CP,1-TP,1-SSCP,1-SXCP,1-HDE,ASSIGNTABLE,NOES,ALLDISTS,NTH,ROWIN,COLIN,2,lv);
   V = or(V,V2);
   ASSIGNTABLE = and(ASSIGNTABLE,V(1:size(ASSIGNTABLE,1),1:size(ASSIGNTABLE,2)));
   nlast = sum(sum(ASSIGNTABLE));
   for(i=1:100)
      [ASSIGNTABLE] = lockdown(MASTER,ASSIGNTABLE,NOES,IALLDISTS, NTH, ROWIN,COLIN);
      [ASSIGNTABLE] = lockdown(MASTER,ASSIGNTABLE,NOES,ALLDISTS, 15, ROWIN,COLIN);
      [MASTER,ASSIGNTABLE]=updateMASTER(MASTER,ASSIGNTABLE,ROWIN,COLIN);[ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN);
      [ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,TP,HDE,ROWIN,COLIN]=reduceAll(ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,TP,HDE,ROWIN,COLIN);[ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN);
      if(sum(sum(ASSIGNTABLE)) == nlast)
         break;
      end
      nlast = sum(sum(ASSIGNTABLE));
   end
end

for(i=[.3 .1])
   if(sum(sum(MASTER))>8)
      break;
   end
   [MASTER,ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,TP,HDE,ROWIN,COLIN]=combine(SXCP,HDE,MASTER,ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,TP,HDE,ROWIN,COLIN,NOES,IALLDISTS,NTH,S1,S2,RDC1,RDC2,RDC1_VECTORS,RDC2_VECTORS,i);
   
   for(q=1:10)
      [ASSIGNTABLE] = lockdown(MASTER,ASSIGNTABLE,NOES,IALLDISTS, NTH,ROWIN,COLIN);
      [ASSIGNTABLE] = lockdown(MASTER,ASSIGNTABLE,NOES,ALLDISTS, 15,ROWIN,COLIN);
      [MASTER,ASSIGNTABLE]=updateMASTER(MASTER,ASSIGNTABLE,ROWIN,COLIN);[ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN);
      [ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,TP,HDE, ROWIN,COLIN]=reduceAll(ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,TP,HDE, ROWIN,COLIN);[ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN);
   end
   S1 = updateTen(MASTER,RDC1,RDC1_VECTORS);
   S2 = updateTen(MASTER,RDC2,RDC2_VECTORS);
   RP1 = NVR_RDC2PROB(ASSIGNTABLE,RDC1,RDC1_VECTORS,S1);
   RP2 = NVR_RDC2PROB(ASSIGNTABLE,RDC2,RDC2_VECTORS,S2);
end



%make the first few assignments
last = sum(sum(ASSIGNTABLE));
while(sum(sum(MASTER))<5)
   if(size(MASTER,1)-sum(sum(MASTER))>0)   
      V = vote5(CP,SXCP,SSCP,TP,HDE,ASSIGNTABLE,NOES,ALLDISTS,NTH,ROWIN,COLIN,1,1);oV = V;
      nV=TP.*CP.*SXCP.*SSCP.*HDE;nV(1:size(ASSIGNTABLE,1),1:size(ASSIGNTABLE,2))=nV(1:size(ASSIGNTABLE,1),1:size(ASSIGNTABLE,2)).*ASSIGNTABLE;
      V(1:size(ASSIGNTABLE,1),1:size(ASSIGNTABLE,2))=V(1:size(ASSIGNTABLE,1),1:size(ASSIGNTABLE,2)).*ASSIGNTABLE;
      [A]=getASS(V,V, ASSIGNTABLE,0,1);A(1:size(ASSIGNTABLE,1),1:size(ASSIGNTABLE,2))=A(1:size(ASSIGNTABLE,1),1:size(ASSIGNTABLE,2)).*ASSIGNTABLE;
      onV=nV;
      V=thresh(oV,7);nV =thresh(onV,10e-8);
      X = and(A(1:size(ASSIGNTABLE,1),1:size(ASSIGNTABLE,2)),and(V(1:size(ASSIGNTABLE,1),1:size(ASSIGNTABLE,2)),nV));
      for(i=1:size(X,1))
         x = find(X(i,:));
         if(length(x)==1)
            ASSIGNTABLE(i,:)=0;
            ASSIGNTABLE(i,x)=1;
         end
      end
   end
   for(q=1:10)
      [ASSIGNTABLE] = lockdown(MASTER,ASSIGNTABLE,NOES,ALLDISTS, 15,ROWIN,COLIN);
      [MASTER,ASSIGNTABLE]=updateMASTER(MASTER,ASSIGNTABLE,ROWIN,COLIN);[ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN);
      [ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,TP,HDE, ROWIN,COLIN]=reduceAll(ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,TP,HDE, ROWIN,COLIN);[ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN);
   end
   if(sum(sum(ASSIGNTABLE))==last)
      break;
   end
   last = sum(sum(ASSIGNTABLE));
end

for(lv = [10 8 6 4 2])
   V = vote5(CP,SXCP,SSCP,TP,HDE,ASSIGNTABLE,NOES,ALLDISTS,NTH,ROWIN,COLIN,2,20);
   V2 = vote5(1-CP,1-TP,1-SSCP,1-SXCP,1-HDE,ASSIGNTABLE,NOES,ALLDISTS,NTH,ROWIN,COLIN,2,20);
   V3 = vote5(CP',SXCP',SSCP',TP',HDE',ASSIGNTABLE,NOES,ALLDISTS,NTH,ROWIN,COLIN,2,20);
   V4 = vote5(1-CP',1-TP',1-SSCP',1-SXCP',1-HDE',ASSIGNTABLE,NOES,ALLDISTS,NTH,ROWIN,COLIN,2,20);
   X = V+V2+V3'+V4';
   X = thresh(X,max(lv,min(nonzeros(X))+1));
   V = and(X,X);
   ASSIGNTABLE = and(ASSIGNTABLE,V(1:size(ASSIGNTABLE,1),1:size(ASSIGNTABLE,2)));
   nlast = sum(sum(ASSIGNTABLE));
   for(i=1:100)
      [ASSIGNTABLE] = lockdown(MASTER,ASSIGNTABLE,NOES,IALLDISTS, NTH, ROWIN,COLIN);
      [ASSIGNTABLE] = lockdown(MASTER,ASSIGNTABLE,NOES,ALLDISTS, 15, ROWIN,COLIN);
      [MASTER,ASSIGNTABLE]=updateMASTER(MASTER,ASSIGNTABLE,ROWIN,COLIN);[ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN);
      [ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,TP,HDE,ROWIN,COLIN]=reduceAll(ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,TP,HDE,ROWIN,COLIN);[ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN);
      if(sum(sum(ASSIGNTABLE)) == nlast)
         break;
      end
      nlast = sum(sum(ASSIGNTABLE));
   end
end


for(i=[.3 .1])
   if(sum(sum(MASTER))>8)
      break;
   end
   [MASTER,ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,TP,HDE,ROWIN,COLIN]=combine(SXCP,HDE,MASTER,ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,TP,HDE,ROWIN,COLIN,NOES,IALLDISTS,NTH,S1,S2,RDC1,RDC2,RDC1_VECTORS,RDC2_VECTORS,i);
   
   for(q=1:10)
      [ASSIGNTABLE] = lockdown(MASTER,ASSIGNTABLE,NOES,IALLDISTS, NTH,ROWIN,COLIN);
      [ASSIGNTABLE] = lockdown(MASTER,ASSIGNTABLE,NOES,ALLDISTS, 15,ROWIN,COLIN);
      [MASTER,ASSIGNTABLE]=updateMASTER(MASTER,ASSIGNTABLE,ROWIN,COLIN);[ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN);
      [ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,TP,HDE, ROWIN,COLIN]=reduceAll(ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,TP,HDE, ROWIN,COLIN);[ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN);
   end
   S1 = updateTen(MASTER,RDC1,RDC1_VECTORS);
   S2 = updateTen(MASTER,RDC2,RDC2_VECTORS);
   RP1 = NVR_RDC2PROB(ASSIGNTABLE,RDC1,RDC1_VECTORS,S1);
   RP2 = NVR_RDC2PROB(ASSIGNTABLE,RDC2,RDC2_VECTORS,S2);
end

[S1 flag1]= updateTen(MASTER,RDC1,RDC1_VECTORS);
if(flag1==0)
   RP1 = NVR_RDC2PROB(ASSIGNTABLE,RDC1,RDC1_VECTORS,S1);
else
   RP1  = (ASSIGNTABLE*0+1);
end

[S2 flag2] = updateTen(MASTER,RDC2,RDC2_VECTORS);
if(flag2==0)
   RP2 = NVR_RDC2PROB(ASSIGNTABLE,RDC2,RDC2_VECTORS,S2);
else
   RP2  = (ASSIGNTABLE*0+1);
end

while(flag1==1 | flag2==1)
   if(size(MASTER,1)-sum(sum(MASTER))>0)   
      V = vote5(CP,SXCP,SSCP,TP,HDE,ASSIGNTABLE,NOES,ALLDISTS,NTH,ROWIN,COLIN,1,1);oV = V;
      nV=TP.*CP.*SXCP.*SSCP.*HDE;nV(1:size(ASSIGNTABLE,1),1:size(ASSIGNTABLE,2))=nV(1:size(ASSIGNTABLE,1),1:size(ASSIGNTABLE,2)).*ASSIGNTABLE;
      
      V(1:size(ASSIGNTABLE,1),1:size(ASSIGNTABLE,2))=V(1:size(ASSIGNTABLE,1),1:size(ASSIGNTABLE,2)).*ASSIGNTABLE;
      [A]=getASS(V,V, ASSIGNTABLE,0,1);A(1:size(ASSIGNTABLE,1),1:size(ASSIGNTABLE,2))=A(1:size(ASSIGNTABLE,1),1:size(ASSIGNTABLE,2)).*ASSIGNTABLE;
      
      onV=nV;
      V=thresh(oV,16);nV =thresh(onV,10e-8);
      X = and(A(1:size(ASSIGNTABLE,1),1:size(ASSIGNTABLE,2)),and(V(1:size(ASSIGNTABLE,1),1:size(ASSIGNTABLE,2)),nV));
      
      for(i=1:size(X,1))
         x = find(X(i,:));
         if(length(x)==1)
            ASSIGNTABLE(i,:)=0;
            ASSIGNTABLE(i,x)=1;
         end
      end
   end
   for(q=1:10)
      [ASSIGNTABLE] = lockdown(MASTER,ASSIGNTABLE,NOES,ALLDISTS, 15,ROWIN,COLIN);
      [MASTER,ASSIGNTABLE]=updateMASTER(MASTER,ASSIGNTABLE,ROWIN,COLIN);[ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN);
      [ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,TP,HDE, ROWIN,COLIN]=reduceAll(ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,TP,HDE, ROWIN,COLIN);[ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN);
   end
   if(sum(sum(ASSIGNTABLE))==last)
      break;
   end
   last = sum(sum(ASSIGNTABLE));
   
   [S1 flag1]= updateTen(MASTER,RDC1,RDC1_VECTORS);
   [S2 flag2] = updateTen(MASTER,RDC2,RDC2_VECTORS);
end


RP1 = NVR_RDC2PROB(ASSIGNTABLE,RDC1,RDC1_VECTORS,S1);
RP2 = NVR_RDC2PROB(ASSIGNTABLE,RDC2,RDC2_VECTORS,S2);

isum = sum(sum(ASSIGNTABLE));
for(lv = 1:10)
   V = vote(CP,SXCP,SSCP,TP,RP1,RP2,HDE,ASSIGNTABLE,NOES,ALLDISTS,NTH,ROWIN,COLIN,2,20);
   V2 = vote(1-CP,1-TP,1-SSCP,1-SXCP,1-RP1,1-RP2,1-HDE,ASSIGNTABLE,NOES,ALLDISTS,NTH,ROWIN,COLIN,2,20);
   V3 = vote(CP',SXCP',SSCP',TP',RP1',RP2',HDE',ASSIGNTABLE,NOES,ALLDISTS,NTH,ROWIN,COLIN,2,20);
   V4 = vote(1-CP',1-TP',1-SSCP',1-SXCP',1-RP1',1-RP2',1-HDE',ASSIGNTABLE,NOES,ALLDISTS,NTH,ROWIN,COLIN,2,20);
   
   X = V+V2+V3'+V4';
   X = thresh(X,max(45,min(nonzeros(X))+40));
   V = and(X,X);
   ASSIGNTABLE = and(ASSIGNTABLE,V(1:size(ASSIGNTABLE,1),1:size(ASSIGNTABLE,2)));
   nlast = sum(sum(ASSIGNTABLE));
   for(i=1:100)
      [ASSIGNTABLE] = lockdown(MASTER,ASSIGNTABLE,NOES,IALLDISTS, NTH, ROWIN,COLIN);
      [ASSIGNTABLE] = lockdown(MASTER,ASSIGNTABLE,NOES,ALLDISTS, 15, ROWIN,COLIN);
      [MASTER,ASSIGNTABLE]=updateMASTER(MASTER,ASSIGNTABLE,ROWIN,COLIN);[ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN);
      [ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,TP,HDE,ROWIN,COLIN]=reduceAll(ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,TP,HDE,ROWIN,COLIN);[ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN);
      if(sum(sum(ASSIGNTABLE)) == nlast)
         break;
      end
      nlast = sum(sum(ASSIGNTABLE));
   end
   if(isum == sum(sum(ASSIGNTABLE)))
      break;
   end
   isum = sum(sum(ASSIGNTABLE));
end


fail=0;
for(i=[.1])
   [MASTER,ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,TP,HDE,ROWIN,COLIN]= ...
       combine(RP1.*RP2,TP,MASTER,ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,TP,HDE,ROWIN,COLIN,NOES,ALLDISTS,15,S1,S2,RDC1,RDC2,RDC1_VECTORS,RDC2_VECTORS, i);
   for(q=1:10)
      [ASSIGNTABLE] = lockdown(MASTER,ASSIGNTABLE,NOES,IALLDISTS, NTH,ROWIN,COLIN);

      [ASSIGNTABLE] = lockdown(MASTER,ASSIGNTABLE,NOES,ALLDISTS, 15,ROWIN,COLIN);
      
      [MASTER,ASSIGNTABLE]=updateMASTER(MASTER,ASSIGNTABLE,ROWIN,COLIN);[ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN);

      [ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,TP,HDE, ROWIN,COLIN]=reduceAll(ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,TP,HDE, ROWIN,COLIN);[ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN);
   
   end
   S1 = updateTen(MASTER,RDC1,RDC1_VECTORS);
   S2 = updateTen(MASTER,RDC2,RDC2_VECTORS);
   RP1 = NVR_RDC2PROB(ASSIGNTABLE,RDC1,RDC1_VECTORS,S1);
   RP2 = NVR_RDC2PROB(ASSIGNTABLE,RDC2,RDC2_VECTORS,S2);
end

currsize=sum(sum(ASSIGNTABLE));
sthr=10;VTH=79;BTH=20;

while(sum(sum(MASTER))<size(MASTER,1))
   if(sum(sum(MASTER))>=size(MASTER,1)-1)
      break;
   end
   nV=TP.*CP.*SXCP.*SSCP.*HDE.*RP1.*RP2;nV(1:size(ASSIGNTABLE,1),1:size(ASSIGNTABLE,2))=nV(1:size(ASSIGNTABLE,1),1:size(ASSIGNTABLE,2)).*ASSIGNTABLE;
   onV=ren(nV);
   
   %do voting
   V = vote(CP,SXCP,SSCP,TP,RP1,RP2,HDE, ASSIGNTABLE,NOES,ALLDISTS,NTH,ROWIN,COLIN,1,lv);oV=V;
   V(1:size(ASSIGNTABLE,1),1:size(ASSIGNTABLE,2))=V(1:size(ASSIGNTABLE,1),1:size(ASSIGNTABLE,2)).*ASSIGNTABLE;
   [A]=getASS(V,V, ASSIGNTABLE,0,1);A(1:size(ASSIGNTABLE,1),1:size(ASSIGNTABLE,2))=A(1:size(ASSIGNTABLE,1),1:size(ASSIGNTABLE,2)).*ASSIGNTABLE;
   
   V=thresh(oV,VTH);nV =thresh(ren(onV),.13);
   X = and(A(1:size(ASSIGNTABLE,1),1:size(ASSIGNTABLE,2)),and(V(1:size(ASSIGNTABLE,1),1:size(ASSIGNTABLE,2)),nV));
   
   while(sum(sum(X))==0 & VTH>0)
      VTH=VTH-1;
      V=thresh(oV,VTH);nV =thresh(ren(onV),.13);
      X = and(A(1:size(ASSIGNTABLE,1),1:size(ASSIGNTABLE,2)),and(V(1:size(ASSIGNTABLE,1),1:size(ASSIGNTABLE,2)),nV));
   end
   
   if(sum(sum(X))==0 & currsize==sum(sum(ASSIGNTABLE)))
      X = thresh(V,max(max(V)));X = X(1:size(ASSIGNTABLE,1),1:size(ASSIGNTABLE,2));X=and(X,X);
   end
   
   for(i=1:size(X,1))
      x = find(X(i,:));
      if(length(x)==1)
         ASSIGNTABLE(i,:)=0;
         ASSIGNTABLE(i,x)=1;
      end
   end
   
   if(sum(sum(MASTER))>=size(MASTER,1)-1)
      break;
   end
   
   for(q=1:10)
      [ASSIGNTABLE] = lockdown(MASTER,ASSIGNTABLE,NOES,ALLDISTS, BTH,ROWIN,COLIN);
      [MASTER,ASSIGNTABLE]=updateMASTER(MASTER,ASSIGNTABLE,ROWIN,COLIN);[ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN);
      [ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,TP,HDE, ROWIN,COLIN]=reduceAll(ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,TP,HDE, ROWIN,COLIN);[ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN);
      S1 = updateTen(MASTER,RDC1,RDC1_VECTORS);
      S2 = updateTen(MASTER,RDC2,RDC2_VECTORS);
      RP1 = NVR_RDC2PROB(ASSIGNTABLE,RDC1,RDC1_VECTORS,S1);
      RP2 = NVR_RDC2PROB(ASSIGNTABLE,RDC2,RDC2_VECTORS,S2);
   end
   
   if(sum(sum(MASTER))>=size(MASTER,1)-1)
      break;
   end
   
   
   for(i=[.4])
      [MASTER,ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,TP,HDE,ROWIN,COLIN]=combine(SXCP,HDE,MASTER,ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,TP,HDE,ROWIN,COLIN,NOES,ALLDISTS,15,S1,S2,RDC1,RDC2,RDC1_VECTORS,RDC2_VECTORS,i);
      for(q=1:10)
         [ASSIGNTABLE] = lockdown(MASTER,ASSIGNTABLE,NOES,IALLDISTS, NTH,ROWIN,COLIN);
         [ASSIGNTABLE] = lockdown(MASTER,ASSIGNTABLE,NOES,ALLDISTS, BTH,ROWIN,COLIN);
         [MASTER,ASSIGNTABLE]=updateMASTER(MASTER,ASSIGNTABLE,ROWIN,COLIN);[ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN);
         [ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,TP,HDE, ROWIN,COLIN]=reduceAll(ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,TP,HDE, ROWIN,COLIN);[ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN);
      end
      S1 = updateTen(MASTER,RDC1,RDC1_VECTORS);
      S2 = updateTen(MASTER,RDC2,RDC2_VECTORS);
      RP1 = NVR_RDC2PROB(ASSIGNTABLE,RDC1,RDC1_VECTORS,S1);
      RP2 = NVR_RDC2PROB(ASSIGNTABLE,RDC2,RDC2_VECTORS,S2);
   end
   
   for(lv = sthr)
      
      V = vote(CP,SXCP,SSCP,TP,RP1,RP2,HDE,ASSIGNTABLE,NOES,ALLDISTS,NTH,ROWIN,COLIN,2,lv);
      V2 = vote(1-CP,1-TP,1-SSCP,1-SXCP,1-RP1,1-RP2,1-HDE,ASSIGNTABLE,NOES,ALLDISTS,NTH,ROWIN,COLIN,2,lv);
      V = or(V,V2);
      ASSIGNTABLE = and(ASSIGNTABLE,V(1:size(ASSIGNTABLE,1),1:size(ASSIGNTABLE,2)));
      nlast = sum(sum(ASSIGNTABLE));
      for(i=1:100)
         [ASSIGNTABLE] = lockdown(MASTER,ASSIGNTABLE,NOES,IALLDISTS, NTH, ROWIN,COLIN);
         [ASSIGNTABLE] = lockdown(MASTER,ASSIGNTABLE,NOES,ALLDISTS, BTH, ROWIN,COLIN);
         [MASTER,ASSIGNTABLE]=updateMASTER(MASTER,ASSIGNTABLE,ROWIN,COLIN);[ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN);
         [ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,TP,HDE,ROWIN,COLIN]=reduceAll(ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,TP,HDE,ROWIN,COLIN);[ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN);
         S1 = updateTen(MASTER,RDC1,RDC1_VECTORS);
         S2 = updateTen(MASTER,RDC2,RDC2_VECTORS);
         RP1 = NVR_RDC2PROB(ASSIGNTABLE,RDC1,RDC1_VECTORS,S1);
         RP2 = NVR_RDC2PROB(ASSIGNTABLE,RDC2,RDC2_VECTORS,S2);
         if(sum(sum(ASSIGNTABLE)) == nlast)
            break;
         end
         nlast = sum(sum(ASSIGNTABLE));
      end
   end
   sthr=max(sthr-1,1);
   
   if(currsize==sum(sum(ASSIGNTABLE)) & VTH>0);
      VTH=VTH-1
   elseif(currsize==sum(sum(ASSIGNTABLE)))
      break
   end
   
   currsize=sum(sum(ASSIGNTABLE));
   BTH=max(NTH,BTH*.9);
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

%HD_SCORE = -1000;
weightedHD_Score = -1000;


HD_SCORE    = computeHD_Score( RDC1, RDC2, RDC1_VECTORS,  ...
			       HDEXCHANGE, HBOND, peakIDs, HSHIFTS, ...
			       NSHIFTS, TYPES, SSTRUCT, NOES, ALLDISTS, IALLDISTS,...
			       SHIFTS_Filename, ...
			       SHIFTX_Filename, MASTER);

%weightedHD_SCORE = computeWeightedHD_Score (RDC1, RDC2, VECTORS, ...
%					    HDEXCHANGE, HBOND, peakIDs, HSHIFTS, ...
%					    NSHIFTS, TYPES, SSTRUCT, NOES, ALLDISTS, IALLDISTS,...
%					    SHIFTS_Filename, SHIFTX_Filename, ...
%					    MASTER);

[assignmentAccuracy,assignments] = computeAssignmentAccuracy(peakIDs, RESNUMS, MASTER);

return









function [ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,TP,HDE,ROWIN,COLIN]=reduceAll(ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,TP,HDE,ROWIN,COLIN);
[ASSIGNTABLE,CP,x,y]=reduce(ASSIGNTABLE,CP,ROWIN,COLIN);
[ASSIGNTABLE,SXCP,x,y]=reduce(ASSIGNTABLE,SXCP,ROWIN,COLIN);
[ASSIGNTABLE,SSCP,x,y]=reduce(ASSIGNTABLE,SSCP,ROWIN,COLIN);
[ASSIGNTABLE,TP,x,y]=reduce(ASSIGNTABLE,TP,ROWIN,COLIN);
[ASSIGNTABLE,HDE,x,y]=reduce(ASSIGNTABLE,HDE,ROWIN,COLIN);
[ASSIGNTABLE,RP1,x,y]=reduce(ASSIGNTABLE,RP1,ROWIN,COLIN);
[ASSIGNTABLE,RP2,x,y]=reduce(ASSIGNTABLE,RP2,ROWIN,COLIN);
[ASSIGNTABLE,ASSIGNTABLE,ROWIN,COLIN]=reduce(ASSIGNTABLE,ASSIGNTABLE,ROWIN,COLIN);
ASSIGNTABLE=and(ASSIGNTABLE,ASSIGNTABLE);


function  [MASTER,ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,TP,HDE, ROWIN,COLIN]=combine(A,B,MASTER,ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,TP,HDE, ROWIN,COLIN,NOES,ALLDISTS,NTH,S1,S2,RDC1,RDC2,RDC1_VECTORS,RDC2_VECTORS,THR);
THR;
if(THR>0)
   [ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN);
   [ASSIGNTABLE] = getUACMB(A,B,ASSIGNTABLE,NOES,ALLDISTS,NTH,THR,1,ROWIN,COLIN,MASTER);[ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN);
   [ASSIGNTABLE] = lockdown(MASTER,ASSIGNTABLE,NOES,ALLDISTS, NTH,ROWIN,COLIN);
   [MASTER,ASSIGNTABLE]=updateMASTER(MASTER,ASSIGNTABLE,ROWIN,COLIN);[ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN);
   [ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,TP,HDE, ROWIN,COLIN]=reduceAll(ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,TP,HDE, ROWIN,COLIN);[ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN);
   S1 = updateTen(MASTER,RDC1,RDC1_VECTORS);
   S2 = updateTen(MASTER,RDC2,RDC2_VECTORS);
   RP1 = NVR_RDC2PROB(ASSIGNTABLE,RDC1,RDC1_VECTORS,S1);
   RP2 = NVR_RDC2PROB(ASSIGNTABLE,RDC2,RDC2_VECTORS,S2);
end



function V = vote(CP,SXCP,SSCP,TP,RP1,RP2,HDE,ASSIGNTABLE,NOES,ALLDISTS,NTH,ROWIN,COLIN,atype,top)

V =CP*0;

if (size(ASSIGNTABLE,1) == 0)
  return
end

for(i=1:7)
   m = nchoosek(1:7,i);
   for(j = 1:size(m,1))
      X = CP*0+1;
      c = m(j,:);
      for(k=1:length(c))
         if(c(k)==1)
            X = X.*CP;
         elseif(c(k)==2)
            X = X.*SXCP;
         elseif(c(k)==3)
            X = X.*SSCP;
         elseif(c(k)==4)
            X = X.*TP;
         elseif(c(k)==5)
	   X = X.*RP1;
         elseif(c(k)==6)
            X = X.*RP2;
         elseif(c(k)==7)
            X = X.*HDE;
         end
      end
      A = getASS(X,X, ASSIGNTABLE,0,atype,top);
      A = A(1:size(V,1),1:size(V,2));
      V = V+A;
   end
end



function V = vote5(CP,SXCP,SSCP,TP,HDE,ASSIGNTABLE,NOES,ALLDISTS,NTH,ROWIN,COLIN,atype,top)
V =CP*0;

if (size(ASSIGNTABLE,1) == 0)
  return
end

for(i=1:5)
   m = nchoosek(1:5,i);
   for(j = 1:size(m,1))
      X = CP*0+1;
      c = m(j,:);
      for(k=1:length(c))
         if(c(k)==1)
            X = X.*CP;
         elseif(c(k)==2)
            X = X.*SXCP;
         elseif(c(k)==3)
            X = X.*SSCP;
         elseif(c(k)==4)
            X = X.*TP;
         elseif(c(k)==5)
            X = X.*HDE;
         end
      end
      A = getASS(X,X, ASSIGNTABLE,0,atype,top);
      A = A(1:size(V,1),1:size(V,2));
      V = V+A;
   end
end


