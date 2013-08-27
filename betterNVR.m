function [assignmentAccuracy,NVR_SCORE,oldNVR_SCORE]= betterNVR(peaks,rdcs,HDEXCHANGE, peakIDs, inputNOES, ...
    inputVECTORS,TYPES, ...
	RESNUMS,SSTRUCT, HBOND, inputALLDISTS, SHIFTS_Filename, SHIFTX_Filename)


%The differences between this and original NVR are:
%This uses also MASTER in lockdown.
% $$$ This does not arbitrarily reduces the bpg's corresponding to SHIFTS, ...
% $$$     SHIFTX, BMRB, according to top m entries.
% $$$ It also does not reduce to a linear number of edges.
% $$$ Takes as input shifts and shiftx filenames. Returns assignment ...
% $$$     accuracy and some NVR scores.
% $$$ This computes NVR score.
% $$$ This does fix some bugs in the original NVR, such as division by 0, ...
% $$$     or restoring too many entries in the matrices if a row is empty. 

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

SHIFTS_Filename
SHIFTX_Filename

HSHIFTS = peaks(:,1);
NSHIFTS = peaks(:,2);


global RDC1 RDC2
global VECTORS finalAssignments
global numPeaks unassignablePeaks
global NTH
global NOES ALLDISTS 
global MASTER

RDC1 = rdcs(:,1);
RDC2 = rdcs(:,2);

VECTORS   = inputVECTORS;

numPeaks          = length(HSHIFTS);
unassignablePeaks = zeros(numPeaks,1);


NOES      = inputNOES;
ALLDISTS  = inputALLDISTS;

%compute a tollerance for NOE distances
NTH=4.8;
mu=mean(mean(ALLDISTS));
if(mu-12.9>0)
   NTH=NTH+(mu-12.9); NTH=min(NTH,8);
end



ASSIGNTABLE = ones(length(HSHIFTS),size(VECTORS,1))/size(VECTORS,1);

%This is the master assignment table
MASTER    = ASSIGNTABLE*0;


%these keep track of which peaks and residues are represented in the current
%matricies

fprintf(1, 'num potential assignments is %d\n',(sum(sum(ASSIGNTABLE))));

ROWIN=1:size(ASSIGNTABLE,1);
COLIN=1:size(ASSIGNTABLE,2);




%prune via amide exchange
HDE = NVR_HD2PROB(ASSIGNTABLE,HDEXCHANGE,HBOND);
ASSIGNTABLE = and(ASSIGNTABLE,HDE);

fprintf(1, 'after HD pruning num potential assignments is %d\n',(sum(sum(ASSIGNTABLE))));


nlast = sum(sum(ASSIGNTABLE));
for(i=1:100)
   [ASSIGNTABLE] = lockdown(MASTER,ASSIGNTABLE,NOES,ALLDISTS, NTH,ROWIN,COLIN);
   if(sum(sum(ASSIGNTABLE )) == nlast)
      break;
   end
   fprintf(1, 'after NOE pruning num potential assignments is %d\n',(sum(sum(ASSIGNTABLE))));
   nlast = sum(sum(ASSIGNTABLE ));
end

%prune via BMRB stats
[CP] = NVR_CS2PROB(ASSIGNTABLE,HSHIFTS,NSHIFTS,TYPES,SSTRUCT,NOES,ALLDISTS,NTH,ROWIN,COLIN);SXCP=CP;SSCP = CP;RP1 = CP;RP2 = CP;
ASSIGNTABLE = and(ASSIGNTABLE,CP);
fprintf(1, 'after CP pruning num potential assignments is %d\n',(sum(sum(ASSIGNTABLE))));

nlast = sum(sum(ASSIGNTABLE));
for(i=1:100)
  [ASSIGNTABLE] = lockdown(MASTER,ASSIGNTABLE,NOES,ALLDISTS, NTH,ROWIN,COLIN);
  if(sum(sum(ASSIGNTABLE)) == nlast)
    break;
  end
  fprintf(1, 'after NOE pruning num potential assignments is %d\n',(sum(sum(ASSIGNTABLE))));
  nlast = sum(sum(ASSIGNTABLE));
end

%prune via the program shiftx
[SXCP] = NVR_SHIFTX2PROB(ASSIGNTABLE,HSHIFTS,NSHIFTS,TYPES,SSTRUCT, ...
			 NOES,ALLDISTS,NTH,ROWIN,COLIN, SHIFTX_Filename);
ASSIGNTABLE = and(ASSIGNTABLE,SXCP);
fprintf(1, 'after SXCP pruning num potential assignments is %d\n',(sum(sum(ASSIGNTABLE))));

   nlast = sum(sum(ASSIGNTABLE));
   for(i=1:100)
      [ASSIGNTABLE] = lockdown(MASTER,ASSIGNTABLE,NOES,ALLDISTS, NTH,ROWIN,COLIN);
      if(sum(sum(ASSIGNTABLE)) == nlast)
         break;
      end
      fprintf(1, 'after NOE pruning num potential assignments is %d\n',(sum(sum(ASSIGNTABLE))));
      
      nlast = sum(sum(ASSIGNTABLE));
   end
% $$$ end



%prune via the program shifts
[SSCP] = NVR_SHIFTS2PROB(ASSIGNTABLE,HSHIFTS,NSHIFTS,TYPES,SSTRUCT, ...
			 NOES,ALLDISTS,NTH,ROWIN,COLIN, SHIFTS_Filename, SHIFTX_Filename);
fprintf(1, 'after SSCP pruning num potential assignments is %d\n',(sum(sum(ASSIGNTABLE))));
ASSIGNTABLE = and(ASSIGNTABLE,SSCP);



[MASTER, ASSIGNTABLE, ROWIN, COLIN, CP, SXCP, SSCP, RP1, RP2, unassignablePeaks] = assignAndNoePruneCycle(MASTER, ...
						  ASSIGNTABLE, ROWIN, COLIN, CP, ...
						  SXCP, SSCP, RP1, RP2,NOES,ALLDISTS, ...
						  NTH, unassignablePeaks);

[ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN, unassignablePeaks);

[MASTER, ASSIGNTABLE, ROWIN, COLIN, CP, SXCP, SSCP, RP1, RP2, unassignablePeaks] = assignAndNoePruneCycle(MASTER, ...
						  ASSIGNTABLE, ROWIN, COLIN, CP, ...
						  SXCP, SSCP, RP1, RP2,NOES,ALLDISTS, ...
						  NTH, unassignablePeaks);

[S1, S2, RP1, RP2] = updateRDCVariables(MASTER, ASSIGNTABLE, RDC1, ...
					RDC2, VECTORS, CP);


[MASTER, ASSIGNTABLE, ROWIN, COLIN, CP, SXCP, SSCP, RP1, RP2, unassignablePeaks] = assignAndNoePruneCycle(MASTER, ...
						  ASSIGNTABLE, ROWIN, COLIN, CP, ...
						  SXCP, SSCP, RP1, RP2,NOES,ALLDISTS, ...
						  NTH, unassignablePeaks);
[numPeaks, numResidues] = size(MASTER);



finalAssignments        = zeros(1000,1000);



[finalAssignments] = recursiveAssign(ASSIGNTABLE, ROWIN, COLIN, MASTER, RP1, RP2, ...
       CP, SXCP, SSCP, RESNUMS, NOES, ALLDISTS, NTH,RDC1,RDC2,VECTORS, ...
		peakIDs, finalAssignments, 0, unassignablePeaks);

keyboard




%function assign(state)
%
%
%  done = 0;
%
%  if (!done)
%    assign the unambiguous peaks.
%    prune.
%    update assigntable, numAssigned, rowin, colin, numAssignedPeaks;
%    if (numAssignedPeaks == numPeaks)
%      done = 1;
%    else
%      relPeakIndex = 1;
%      numUnassignedPeaks = totalNumPeaks -  numAssignedPeaks;
%
%      while (relPeakIndex < numUnassignedPeaks)
%        assign peak # relPeakIndex;
%        set the new state parameters, ASSIGNTABLE, ROWIN,COLIN, bpgs...
%        ASSIGN(STATE);
%      end while
%
%    end if
%  end if
%  print out the assignment;
%  return;


function [finalAssignments] = recursiveAssign(ASSIGNTABLE, ROWIN, COLIN, MASTER, RP1, RP2, ...
						  CP, SXCP, SSCP, RESNUMS, NOES, ALLDISTS, NTH, ...
						  RDC1,RDC2,VECTORS, ...
						  peakIDs, finalAssignments, ...
						  depth, unassignablePeaks)

done = 0;
depth = depth + 1;

numAssignedPeaks                  = sum(sum(MASTER));
[totalNumPeaks, dummyVar]         = size(MASTER);

fprintf('at depth %d, numAssignedPeaks = %d\n', depth, numAssignedPeaks);

if (numAssignedPeaks == totalNumPeaks)
  [finalAssignments] = computeCorrectness(MASTER, peakIDs, RESNUMS, finalAssignments);
else

  [MASTER, CP, SSCP, SXCP, RP1, RP2, ASSIGNTABLE, numUnassignedPeaks, ROWIN, COLIN, unassignablePeaks] = ...
      doTheUnambiguousAssignment(MASTER, CP, SSCP, SXCP, RP1, RP2, ...
				 ASSIGNTABLE, RDC1, RDC2,VECTORS, ...
				 ROWIN, COLIN, NOES, ALLDISTS, NTH, ...
				 unassignablePeaks);
  
  doTheBranching(MASTER, CP, SXCP, SSCP, RP1, RP2, ROWIN, ...
		 COLIN, numUnassignedPeaks, unassignablePeaks, ...
		 ASSIGNTABLE, RESNUMS);

end %if

return

function [ASSIGNTABLE] = lockdown(MASTER, ASSIGNTABLE,NOES,ALLDISTS, NTH,ROWIN,COLIN)
%prune via NOEs
if (1)
  myROWIN             = 1:size(MASTER,1);
  myCOLIN             = 1:size(MASTER,2);
  MASTER(ROWIN,COLIN) = ASSIGNTABLE;
  NP                  = NVR_NOE2PROB(MASTER,NOES,ALLDISTS,NTH, myROWIN,myCOLIN);
  NP                  = NP(ROWIN,COLIN);
else
  NP = NVR_NOE2PROB(ASSIGNTABLE,NOES,ALLDISTS,NTH,ROWIN,COLIN);
end
%MSA-SB

ASSIGNTABLE           = and(NP,NP);
fprintf(1, 'after NOE pruning num potential assignments is %d\n',(sum(sum(ASSIGNTABLE))));
[sizeAssignTableX, sizeAssignTableY] = size(ASSIGNTABLE);
numCombinations       = 1;

for i = 1:sizeAssignTableX
  numCombinations = numCombinations * length(find(ASSIGNTABLE(i,:)));
end

fprintf(1, 'after NOE pruning num potential combinations is %d\n',numCombinations);


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
if((size(M,1)==5) & (size(M,2)==5) & (rank(M) == 5))
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





function [MASTER,A, terminate, unassignablePeaks]= ...
    updateMASTER(MASTER,A,ROWIN,COLIN, unassignablePeaks);
%first, propagate contraints
recur=1;
terminate = 0;
rm=ones(1,size(A,1));
cm=ones(1,size(A,2));
ain=A;
while(recur==1)
   recur=0;
   for(i=1:size(A,1))
      x=find(A(i,:));
      if(length(x)==1)
	%SB: should we do this? should we erase all others? lets
        %keep for now.
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
      if (length(x)==1 & size(A,1)==size(A,2))
         A(x,:)=0;
         A(x,i)=1;
         if(cm(i)==1)
            cm(i)=0;
            recur=1;
         end
      end
   end
end

% SB should we keep this?;

for(i=1:size(A,1))
   if(sum(A(i,:))==0)
     %let's keep this also, for now.
     %      A(i,:)=ain(i,:);
     %       terminate = 1; %no restoration
% $$$       x=find(A(i,:));
% $$$       A(x,:)=ain(x,:);
     fprintf('could not find a residue for peak %d\n', ROWIN(i));
%     MASTER (ROWIN(i), :) = 1; %assigning it to all 1's, denoting
%                               %that we cannot assign this peak.
     unassignablePeaks(ROWIN(i)) = 1;
     keyboard
     %       return;
   else
     x = find(A(i,:));
     if(length(x)==1)
% $$$      if (ROWIN(i) ~= COLIN(x)) %MSA-DEBUG
% $$$        fprintf(1, 'making a wrong assignment\n');
% $$$        keyboard
% $$$      end
       fprintf('assigning %d to %d in updateMASTER\n', ROWIN(i),COLIN(x));  
       MASTER(ROWIN(i),COLIN(x))=1;
       if (ROWIN(i) ~= COLIN(x))
	 fprintf(1, 'wrong assignment\n');
	 %       keyboard
       end
     end
   end
end

% need also to take care of the peak that is unassigbanle and remove ...
%     thaty
% from OROWIN, COLIN, CP, etc.?


function [ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,ROWIN,COLIN]= ...
    reduceAll(ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,ROWIN,COLIN, unassignablePeaks);
[ASSIGNTABLE,CP,x,y]=reduce(ASSIGNTABLE,CP,ROWIN,COLIN, unassignablePeaks);
[ASSIGNTABLE,SXCP,x,y]=reduce(ASSIGNTABLE,SXCP,ROWIN,COLIN, unassignablePeaks);
[ASSIGNTABLE,SSCP,x,y]=reduce(ASSIGNTABLE,SSCP,ROWIN,COLIN, unassignablePeaks);
[ASSIGNTABLE,ASSIGNTABLE,ROWIN,COLIN]=reduce(ASSIGNTABLE,ASSIGNTABLE,ROWIN,COLIN, unassignablePeaks);
ASSIGNTABLE=and(ASSIGNTABLE,ASSIGNTABLE);

function [A,B,ROWIN,COLIN]=reduce(A,B,ROWIN,COLIN, unassignablePeaks)
ain=A;
bin=B;
ct=1;row=[1:size(A,1)];col=[1:size(A,2)];
for (i=1:size(A,1))
   x=find(A(i,:));
   if (length(x)==1)
      row(i)=0;
      col(x)=0;   
   elseif (length(x) == 0)
       row(i) = 0;
   end
end
B=B.*A;

listOfUnassignablePeaks = find(unassignablePeaks);
if (~isempty(listOfUnassignablePeaks))
  relIndicesOfUnassignablePeaks = find(ROWIN == listOfUnassignablePeaks)
  row(relIndicesOfUnassignablePeaks) = 0;
end

% $$$ for(i=1:size(B,1))
% $$$    if(sum(B(i,:))==0)
% $$$      %      B(i,:)=bin(i,:);
% $$$      row(i) = 0;
% $$$      %   end
% $$$    else
% $$$      %if(sum(B(i,:))==0)
% $$$      % B(i,:)=1;
% $$$      %end
% $$$      B(i,:)=B(i,:)/sum(B(i,:));
% $$$    end
% $$$ end
if(isempty(row))
   return
end
row=find(row);
col=find(col);
B=B(row,col);
for(i=1:size(B,1))
   if (sum(B(i,:))==0)
%SB: this part of the code needs to be checked.
     %     keyboard
     if (size(B,2) == size(bin,2))
       B(i,:)=bin(i,:);
     else
       B(i,:) = 1;
     end
   end
   B(i,:)=B(i,:)/sum(B(i,:));
end
ROWIN=ROWIN(row);
COLIN=COLIN(col);


%pulls unassigned peak and residue indices to the front. bu uses a
%crude way of checking whether a peak or residu is assigned. maybe
%there is a peak that is unassignable? what to do in theat case?

function [ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN, unassignablePeaks)
sct=1;

sizeROWINATENTRANCE = length(ROWIN);

for(i=1:size(MASTER,1))
   x = find(MASTER(i,:));
   if ((length(x)==0) & (unassignablePeaks(i) == 0)) %corresponds to an unassigned, but assignable, peak.
      ROWIN(sct)=i;
%SB
      %only if the peak has been assigned rowin is reduces. halbuki it must ...
%      be reduced also if that peak is unassignable as determined in ...
%      assigntable. Maybe we can assign to everything in MASTER to ...
%      denote that.
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

SIZErowinatExit = length(ROWIN);
if (SIZErowinatExit ~= sizeROWINATENTRANCE) 
  fprintf('here the problem lise.');
  keyboard
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
     sum_ithRow = sum(HM(i,:));
     if (sum_ithRow == 0)
       HM(i,:) = 1/size(HM,2);
     else
       HM(i,:)=HM(i,:)/sum_ithRow;
     end
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

		      
function M =ren(M)

for(i=1:size(M))
   if(sum(M(i,:))==0)
      M(i,:)=1;
   end
   M(i,:)=M(i,:)/sum(M(i,:));   
end


function [V,numVoters] = getVotes(CP,SXCP,SSCP,ASSIGNTABLE,numAssignedPeaks, MASTER,RDC1,RDC2,VECTORS)
V1 = getASS(CP  ,CP  ,ASSIGNTABLE, 0, 1, 0);
V2 = getASS(SSCP,SSCP,ASSIGNTABLE, 0, 1, 0);
V3 = getASS(SXCP,SXCP,ASSIGNTABLE, 0, 1, 0);
V  = V1 + V2 + V3;
numVoters = 3;

if (numAssignedPeaks >= 5)
  S1  = ones(3);S2=S1;
  S1  = updateTen(MASTER,S1,RDC1,VECTORS);
  S2  = updateTen(MASTER,S2,RDC2,VECTORS);
  RP1 = NVR_RDC2PROB(ASSIGNTABLE,RDC1,VECTORS,S1);
  RP2 = NVR_RDC2PROB(ASSIGNTABLE,RDC2,VECTORS,S2);
  V4  = getASS      (RP1, RP1, ASSIGNTABLE, 0, 1, 0);
  V5  = getASS      (RP2, RP2, ASSIGNTABLE, 0, 1, 0);
  V   = V + V4 + V5;
  numVoters = 5;
end

V = thresh(V, numVoters);
[numRemainingPeaks, numRemainingResidues] = size(ASSIGNTABLE);
V = V(1:numRemainingPeaks, 1:numRemainingResidues);

function [bestAssignmentPeak, bestAssignmentResidue] = pickBestAssignmentPeakAndResidue(V, CP, SXCP, SSCP, RP1, RP2, numVoters);    

maxProbabilityOfAssignment = 0;
bestAssignmentPeak = -1;
bestAssignmentResidue = -1;
for (i=1:size(V,1))
  
  x = find(V(i,:));
  
  if (length(x)==1)
    probabilityOfAssignment = CP(i,x) * SXCP (i,x) * SSCP(i,x);
    if (numVoters == 5)
      probabilityOfAssignment = probabilityOfAssignment * RP1(i,x) * RP2(i,x);
    end
    
    if (probabilityOfAssignment > maxProbabilityOfAssignment)
      maxProbabilityOfAssignment = probabilityOfAssignment;
      bestAssignmentPeak    = i;
      bestAssignmentResidue = x;
    end
    %	ASSIGNTABLE(i,:)=0;
    %	ASSIGNTABLE(i,x)=1;
  end
end

if (bestAssignmentResidue ~= -1)
  fprintf(1, 'found th best assignment shuld be done with probability %f\n', maxProbabilityOfAssignment);
end


function [MASTER, ASSIGNTABLE, ROWIN, COLIN, CP, SXCP, SSCP, RP1, RP2, unassignablePeaks] = assignAndUpdateBPGs(MASTER, ASSIGNTABLE, ROWIN, COLIN, CP, SXCP, SSCP, RP1, RP2, unassignablePeaks);


[MASTER,ASSIGNTABLE,terminate, unassignablePeaks]=updateMASTER(MASTER,ASSIGNTABLE, ...
						  ROWIN,COLIN, unassignablePeaks);
[ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN, unassignablePeaks);
if (terminate)
  return;
end

[ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,ROWIN,COLIN]= ...
    reduceAll(ASSIGNTABLE,CP,SXCP,SSCP,RP1,RP2,ROWIN,COLIN, unassignablePeaks);
[ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN, unassignablePeaks);


function [S1, S2, RP1, RP2] = updateRDCVariables(MASTER, ASSIGNTABLE, RDC1, ...
						 RDC2, VECTORS, CP)

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


function [sortedDistinction, sortedRelPeakIndices] = findPeakToBranchOn(MASTER, numUnassignedPeaks, CP, SSCP, SXCP, ...
		       RP1, RP2);

distinction = zeros(numUnassignedPeaks,1);

for i = 1:numUnassignedPeaks
  weights = CP (i,:) .* ...
      SXCP (i, :) .* SSCP ...
      (i, :);
  if (sum(sum(MASTER)) >= 5)
    weights = weights .* RP1 (i, :) .* ...
	RP2(i, :);
  end
  
  sortedWeights = sort(-weights);
  distinction(i)  = abs(sortedWeights(1) - ...
			sortedWeights(2));
end

[sortedDistinction, sortedRelPeakIndices] = sort(-distinction);
sortedDistinction = -sortedDistinction;

function [sortedResidueProbabilities, relResidueIndices, numCandidateResidues] = findTheOrderInWhichResiduesAreGoingToBeTried(CP, SXCP, SSCP, RP1, RP2, relPeakIndex, ...
						  ASSIGNTABLE, MASTER);


      
residueCandidatesRelIndices = find(ASSIGNTABLE(relPeakIndex,:));
numCandidateResidues = length(residueCandidatesRelIndices);

AGGREGATE_VOTE = CP (relPeakIndex,residueCandidatesRelIndices) .* ...
    SXCP (relPeakIndex, residueCandidatesRelIndices) .* SSCP ...
    (relPeakIndex, residueCandidatesRelIndices);
if (sum(sum(MASTER)) >= 5)
  AGGREGATE_VOTE = AGGREGATE_VOTE .* RP1 (relPeakIndex, residueCandidatesRelIndices) .* ...
      RP2(relPeakIndex, residueCandidatesRelIndices);
end


[sortedResidueProbabilities,sortedRelResidueIndicesWRT_residueCandidatesRelIndices] ...
    = sort(-AGGREGATE_VOTE);

residueCandidatesRelIndices;
AGGREGATE_VOTE;
sortedResidueProbabilities
%	  keyboard
relResidueIndices = residueCandidatesRelIndices(sortedRelResidueIndicesWRT_residueCandidatesRelIndices);


function [MASTER, ASSIGNTABLE, ROWIN, COLIN, CP, SXCP, SSCP, RP1, ...
	  RP2, unassignablePeaks] = assignAndNoePruneCycle(MASTER, ASSIGNTABLE, ROWIN, COLIN, CP, SXCP, SSCP, RP1, RP2,NOES,ALLDISTS, NTH, unassignablePeaks);

nlast = sum(sum(ASSIGNTABLE));

for(i=1:100)
  
  [MASTER, ASSIGNTABLE, ROWIN, COLIN, CP, SXCP, SSCP, RP1, RP2, unassignablePeaks] = assignAndUpdateBPGs(MASTER, ...
						  ASSIGNTABLE, ROWIN, COLIN, CP, ...
						  SXCP, SSCP, RP1, RP2,unassignablePeaks);

  [ASSIGNTABLE] = lockdown(MASTER,ASSIGNTABLE,NOES,ALLDISTS, ...
			   NTH,ROWIN,COLIN);

  fprintf(1, 'after NOE pruning num potential assignments is %d\n',(sum(sum(ASSIGNTABLE))));
  
  [MASTER, ASSIGNTABLE, ROWIN, COLIN, CP, SXCP, SSCP, RP1, RP2, unassignablePeaks] = assignAndUpdateBPGs(MASTER, ...
						  ASSIGNTABLE, ROWIN, COLIN, CP, ...
						  SXCP, SSCP, RP1, ...
						  RP2, unassignablePeaks);
  if(sum(sum(ASSIGNTABLE)) == nlast)
    break;
  end
  nlast = sum(sum(ASSIGNTABLE));
end


function [finalAssignments] = computeCorrectness(MASTER, peakIDs, RESNUMS, finalAssignments)

numAssignments = numAssignments + 1;
a = load('answerkey.txt');
assignments=0;
correct=0;
incorr=0;

for(i=1:size(MASTER,1))
  pk=peakIDs(i);%get the peak id
  x = find(MASTER(i,:));
  if (isempty(x))
    x = -1;
    %     rn = -1;
    rn = 1000;
  else
    rn=RESNUMS(x);%get the residue it was assigned to 
  end
  
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
assignmentAccuracy = correct/size(MASTER,1)*100;
fprintf('Assignment Accuracy = %f percent \n',assignmentAccuracy);
fprintf('\n');
fprintf('\n');
fprintf('Assignments\n');
fprintf('Peak ID -> Residue Number\n');
for(i=1:size(assignments,1))
  
  finalAssignments (assignments(i,1),assignments(i,2)) = finalAssignments(assignments(i,1),assignments(i,2)) + 1;      
  if(incorr(i)==0)
    fprintf('%d	%d\n',assignments(i,1),assignments(i,2));
  else
    fprintf('*%d	%d  (correct=%d->%d)\n',assignments(i,1), ...
	    assignments(i,2),a(i,1),a(i,2));
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


[finalAssignmentsSize1, finalAssignmentsSize2] = size(finalAssignments);

for i = 1:finalAssignmentsSize1
  residueList = find(finalAssignments(i,:));
  if (~isempty(residueList))
    [i residueList  finalAssignments(i, residueList)*1.0/numAssignments]
  end
end
keyboard

function [MASTER, CP, SSCP, SXCP, RP1, RP2, ASSIGNTABLE, numUnassignedPeaks, ROWIN, COLIN, unassignablePeaks, numAssignedPeaks] = doTheUnambiguousAssignment(MASTER, CP, SSCP, SXCP, RP1, RP2, ASSIGNTABLE, RDC1, RDC2,VECTORS, ROWIN, COLIN, NOES, ALLDISTS, NTH, unassignablePeaks);


[totalNumPeaks, dummyVar] = size(MASTER);
numAssignedPeaks = sum(sum(MASTER));
[V,numVoters] = getVotes(CP,SXCP,SSCP,ASSIGNTABLE,numAssignedPeaks, MASTER,RDC1,RDC2,VECTORS);
  
[bestAssignmentPeak, bestAssignmentResidue] = pickBestAssignmentPeakAndResidue(V, CP, SXCP, SSCP, RP1, RP2, numVoters);    
  
if (bestAssignmentResidue ~= -1)
  ASSIGNTABLE(bestAssignmentPeak,:) = 0;
  ASSIGNTABLE(bestAssignmentPeak, bestAssignmentResidue) = 1;
  
  fprintf(1, 'making the unambigous assignment of %d to %d\n', ROWIN(bestAssignmentPeak), COLIN(bestAssignmentResidue));
end

[MASTER, ASSIGNTABLE, ROWIN, COLIN, CP, SXCP, SSCP, RP1, RP2, unassignablePeaks] = assignAndNoePruneCycle(MASTER, ...
						  ASSIGNTABLE, ROWIN, COLIN, CP, ...
						  SXCP, SSCP, RP1, RP2,NOES,ALLDISTS, ...
						  NTH, unassignablePeaks);

[ROWIN,COLIN]=sanitycheck(MASTER,ROWIN,COLIN, unassignablePeaks); 
%why is this here
%in the middle?
								  
								  
[MASTER, ASSIGNTABLE, ROWIN, COLIN, CP, SXCP, SSCP, RP1, RP2, unassignablePeaks] = assignAndNoePruneCycle(MASTER, ...
						  ASSIGNTABLE, ROWIN, COLIN, CP, ...
						  SXCP, SSCP, RP1, RP2,NOES,ALLDISTS, ...
						  NTH, unassignablePeaks);

[S1, S2, RP1, RP2] = updateRDCVariables(MASTER, ASSIGNTABLE, RDC1, ...
					RDC2, VECTORS);

numAssignedPeaks = sum(sum(MASTER));

fprintf('numAssignedPeaks = %d\n', numAssignedPeaks);
if (numAssignedPeaks == totalNumPeaks)
  done = 1;
else
  numUnassignedPeaks = totalNumPeaks -  numAssignedPeaks;
  %numCandidatesPerPeak = zeros(numUnassignedPeaks, 1);
  
  [debugX,debugY] = size(ASSIGNTABLE);
  if (debugX ~= numUnassignedPeaks) 
    keyboard
  end
end

  
function doTheBranching(MASTER, CP, SXCP, SSCP, RP1, RP2, ROWIN, ...
			COLIN, numUnassignedPeaks, unassignablePeaks, ...
			ASSIGNTABLE, RESNUMS, NOES)

[sortedDistinction, sortedRelPeakIndices] = ...
    findPeakToBranchOn(MASTER, numUnassignedPeaks, CP, SSCP, SXCP, ...
		       RP1, RP2);
numTriedPeaks      =  1;

%      while (relPeakIndex <= numUnassignedPeaks)
while(numTriedPeaks <= numUnassignedPeaks)

  %assign peak # relPeakIndex;
  
  %sorted peaks in terms of their distinction and try them
  %one by one using numTriedPeaksIndex.
  
  %[numCandidates,relPeakIndex] = min(numCandidatesPerPeak);
  %VERY_BIG_NUMBER = 10000;
  %	fprintf(1, 'will try to assign %d whichhas %d
  %	candidates\n',ROWIN(relPeakIndex),numCandidates);
  
  relPeakIndex = sortedRelPeakIndices(numTriedPeaks);
  
  fprintf(1, 'will try to assign %d whichhas %f distinction\n',ROWIN(relPeakIndex),sortedDistinction(numTriedPeaks));
					     
  %numCandidatesPerPeak(relPeakIndex) = VERY_BIG_NUMBER;
					     
  [sortedResidueProbabilities, relResidueIndices, numCandidateResidues] = ...
      findTheOrderInWhichResiduesAreGoingToBeTried(CP, SXCP, SSCP, RP1, RP2,...
						   relPeakIndex, ...
						   ASSIGNTABLE, MASTER);
  COLIN(relResidueIndices)
  
  for candidateResidueIndex = 1:numCandidateResidues
    
    newASSIGNTABLE    = ASSIGNTABLE;
    
    newASSIGNTABLE(relPeakIndex, :) = 0; 
    %	    newASSIGNTABLE(relPeakIndex,
    %	    residueCandidatesRelIndices(relResidueIndex)) ...
    relResidueIndex                    = relResidueIndices(candidateResidueIndex);
    
    newASSIGNTABLE(relPeakIndex, relResidueIndex)  = 1;
    fprintf('branching by assigning %d to %d. prob is %f\n', ROWIN(relPeakIndex), COLIN(relResidueIndex),sortedResidueProbabilities(candidateResidueIndex) );
    
    [sizeNewAssignTableX, sizeNewAssignTableY] = size(newASSIGNTABLE);
    for sizeNewAssignTableX_Counter = 1:sizeNewAssignTableX
      if (isempty(find(newASSIGNTABLE(sizeNewAssignTableX_Counter,:))))
	fprintf('the %d peak has already no possible residue',sizeNewAssignTableX_Counter);
	keyboard
      end
    end
    
    [finalAssignments] = recursiveAssign(newASSIGNTABLE, ROWIN, COLIN, MASTER,  RP1, RP2, ...
						  CP, SXCP, SSCP, RESNUMS, NOES, ALLDISTS, NTH, ...
						  RDC1,RDC2,VECTORS, ...
						  peakIDs, ...
						  finalAssignments, ...
						  depth, unassignablePeaks, RESNUMS);
  end
  relPeakIndex = relPeakIndex + 1;
  numTriedPeaks = numTriedPeaks + 1;
end % while
