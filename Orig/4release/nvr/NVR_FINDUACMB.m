function M = NVR_FINDUACMB(A,B,ASSIGNTABLE,THR)

%NVR_FINDUACMB: This computes assignment probabilities 

%////////////////////////////////////////////////////////////////////////////////////////////
%//  NVR_FINDUACMB.m
%//
%//  Version:		0.1
%//
%//  Description:	 This computes assignment probabilities
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

%    NVR_FINDUACMB
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

%Take in the 3 assignment matricies based on chem shift, return
%as unambiguous, anything that all 3 agree on. 
if(size(ASSIGNTABLE,1)-size(ASSIGNTABLE,2)~=0)
   F=ones(max(size(ASSIGNTABLE)));
   F(1:size(ASSIGNTABLE,1),1:size(ASSIGNTABLE,2))=ASSIGNTABLE;
   ASSIGNTABLE=F;   
   for(i=1:size(ASSIGNTABLE,1))
      ASSIGNTABLE(i,:)=ASSIGNTABLE(i,:)/sum(ASSIGNTABLE(i,:));
   end
end
ASSIGNTABLE=and(ASSIGNTABLE,ASSIGNTABLE);

D = A.*B;
for(i=1:size(D,1))
   D(i,:)=D(i,:)/sum(D(i,:));
end
%take all 3
M = getASS(A,ASSIGNTABLE);
M=and(M,getASS(B,ASSIGNTABLE));


M=M(1:size(A,1),1:size(A,2)).*D;
M=thresh(M,THR);
M=uthresh(M,1);
M=and(M,M);

[length(find(A==1)) sum(diag(M)) sum(sum(M))] ;


return



function M = getASS(X,ASSIGNTABLE) 

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
   h=hungarian(HM*-1);
%   h=simp(HM);

for(i=1:length(h))
   M(i,h(i)) = 1 * ASSIGNTABLE(i,h(i));
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
