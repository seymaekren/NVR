function [M] = HD_CS2PROB(TABLE,H,N,TYPES,SSTRUCT,NOES,ALLDISTS,NTH,ROWIN,COLIN);

%NVR_CS2PROB: This computes assignment probabilities based on BMRB stats, it is not meant
%             to be called by the user


%////////////////////////////////////////////////////////////////////////////////////////////
%//  HD_CS2PROB.m
%//
%//  Version:		0.1
%//
%//  Description:	 This computes assignment probabilities based on BMRB stats
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

%    HD_CS2PROB
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

FILTERS=load('~/NVR/CHEMSHIFTSTATS/FILTERS');FILTERS=FILTERS.FILTERS;

%the stats for the shifts were obtained from the wishart group webpage
%http://redpoll.pharmacy.ualberta.ca/RefDB/stat.html
M = TABLE*0;
[rows,cols]=find(TABLE==1);
if(length(rows)==0)
   rows=-1;
   cols=-1;
end
mH=mean(H);sH=std(H);
mN=mean(N);sN=std(N);

for(i=1:size(TABLE,1))
   h = H(i);
   n = N(i);
   p1 = normpdf(h,mH,sH);
   p2 = normpdf(n,mN,sN);
   prior=p1*p2;
   for(j=1:length(TYPES))
      if(strcmp(TYPES(j),'ALA')==1)
         M(i,j)=getProb(h,n, SSTRUCT(j),FILTERS,1)/prior;
      elseif(strcmp(TYPES(j),'CYS')==1)
         M(i,j)=getProb(h,n,SSTRUCT(j),FILTERS,2)/prior;
      elseif(strcmp(TYPES(j),'ASP')==1)
         M(i,j)=getProb(h,n,SSTRUCT(j),FILTERS,3)/prior;
      elseif(strcmp(TYPES(j),'GLU')==1)
         M(i,j)=getProb(h,n,SSTRUCT(j),FILTERS,4)/prior;
      elseif(strcmp(TYPES(j),'PHE')==1)
         M(i,j)=getProb(h,n,SSTRUCT(j),FILTERS,5)/prior;
      elseif(strcmp(TYPES(j),'GLY')==1)
         M(i,j)=getProb(h,n,SSTRUCT(j),FILTERS,6)/prior;
      elseif(strcmp(TYPES(j),'HIS')==1)
         M(i,j)=getProb(h,n,SSTRUCT(j),FILTERS,7)/prior;
      elseif(strcmp(TYPES(j),'ILE')==1)
         M(i,j)=getProb(h,n,SSTRUCT(j),FILTERS,8)/prior;
      elseif(strcmp(TYPES(j),'LYS')==1)
         M(i,j)=getProb(h,n,SSTRUCT(j),FILTERS,9)/prior;
      elseif(strcmp(TYPES(j),'LEU')==1)
         M(i,j)=getProb(h,n,SSTRUCT(j),FILTERS,10)/prior;
      elseif(strcmp(TYPES(j),'MET')==1)
         M(i,j)=getProb(h,n,SSTRUCT(j),FILTERS,11)/prior;
      elseif(strcmp(TYPES(j),'ASN')==1)
         M(i,j)=getProb(h,n,SSTRUCT(j),FILTERS,12)/prior;%
      elseif(strcmp(TYPES(j),'GLN')==1)
         M(i,j)=getProb(h,n,SSTRUCT(j),FILTERS,14)/prior;
      elseif(strcmp(TYPES(j),'ARG')==1)
         M(i,j)=getProb(h,n,SSTRUCT(j),FILTERS,15)/prior;
      elseif(strcmp(TYPES(j),'SER')==1)
         M(i,j)=getProb(h,n,SSTRUCT(j),FILTERS,16)/prior;
      elseif(strcmp(TYPES(j),'THR')==1)
         M(i,j)=getProb(h,n,SSTRUCT(j),FILTERS,17)/prior;
      elseif(strcmp(TYPES(j),'VAL')==1)
         M(i,j)=getProb(h,n,SSTRUCT(j),FILTERS,18)/prior;
      elseif(strcmp(TYPES(j),'TRP')==1)
         M(i,j)=getProb(h,n,SSTRUCT(j),FILTERS,19)/prior;
      elseif(strcmp(TYPES(j),'TYR')==1)
         M(i,j)=getProb(h,n,SSTRUCT(j),FILTERS,20)/prior;
      else
         PROBLEM = TYPES(j)   
      end
      if ((M(i,j) == 0) & (i == j))
	fprintf(1, 'i = %d j = %d h = %f n = %f\n',i,j,h,n);
	keyboard
      end
   end
   M(i,:) = M(i,:)/sum(M(i,:));%re-normalize
end
M = M .* TABLE;%for HD TABLE is uniform therefore I am not sure
               %whether that makes any difference.

%renornmalize
for(i=1:size(M,1))
   M(i,:)=M(i,:)/sum(M(i,:));
end

for(q=1:10)
   
  M_beforeThresholding                    = M;
  
  M                                       = thresh(M,.01*median(nonzeros(M)));
   
  [entriesSetToZero_X,entriesSetToZero_Y] = find(M- M_beforeThresholding);
  
  fprintf(1, '%d entries have been zeroed in the BMRB matrix.\n', length(entriesSetToZero_X));
  
  M_beforeNOE_Applying = M;
  
  TABLE = and(M,M);
   nlast = sum(sum(TABLE));
   for(i=1:100)
      NP = NVR_NOE2PROB(TABLE(1:size(M,1),:),NOES,ALLDISTS,NTH,ROWIN,COLIN);
      %note noe pruning called here.
      TABLE(1:size(M,1),:)=and(TABLE(1:size(M,1),:),NP);
      if(sum(sum(TABLE)) == nlast)
         break;
      end
      nlast = sum(sum(TABLE));
   end
   
   M = M.*TABLE;

   [entriesNOE_PrunedX, entriesNOE_PrunedY] = find(M-M_beforeNOE_Applying);
   
   fprintf(1, '%d entries have been NOE pruned in the BMRB matrix.\n', length(entriesNOE_PrunedX));
   
%   keyboard
   
   for(i=1:size(M,1))
      M(i,:)=M(i,:)/sum(M(i,:));
   end
end

%renornmalize
for(i=1:size(M,1))
   if (M(i,i) == 0)
     fprintf(1, 'i = %d M(i,i) = 0\n',i);
     keyboard
   end
   M(i,:)=M(i,:)/sum(M(i,:));
end

function p = getProb(h,n,SSTYPE,FILTERS,TY)%


if(strcmp(SSTYPE,'C')==1)
   
   xH = FILTERS(3,TY,1,1);
   xHS = FILTERS(3,TY,1,2);
   xN = FILTERS(3,TY,2,1);
   xNS = FILTERS(3,TY,2,2);
   xHMX = FILTERS(3,TY,1,4)*30;
   xHMN = FILTERS(3,TY,1,3)*30;
   xNMX = FILTERS(3,TY,2,4)*50;%
   xNMN = FILTERS(3,TY,2,3)*50;%
   
elseif(strcmp(SSTYPE,'B')==1)
   
   xH = FILTERS(2,TY,1,1);
   xHS = FILTERS(2,TY,1,2);
   xN = FILTERS(2,TY,2,1);
   xNS = FILTERS(2,TY,2,2);
   xHMX = FILTERS(2,TY,1,4)*13;%
   xHMN = FILTERS(2,TY,1,3)*13;%
   xNMX = FILTERS(2,TY,2,4)*80;%
   xNMN = FILTERS(2,TY,2,3)*80;%
   
else
   xH = FILTERS(1,TY,1,1);
   xHS = FILTERS(1,TY,1,2);
   xN = FILTERS(1,TY,2,1);
   xNS = FILTERS(1,TY,2,2);
   xHMX = FILTERS(1,TY,1,4)*9;
   xHMN = FILTERS(1,TY,1,3)*9;
   xNMX = FILTERS(1,TY,2,4)*30;%
   xNMN = FILTERS(1,TY,2,3)*30;%
   
end

if(h>xH)
   if((h-xH)/xHS > xHMX)
      p=0;return
   else
      p1=normpdf(h,xH,xHS);
   end
else
   if((xH-h)/xHS > xHMN)
      p=0;return
   else
      p1=normpdf(h,xH,xHS);
   end
end
if(n>xN)
   if((n-xN)/xNS > xNMX)
      p=0;return
   else
      p2=normpdf(n,xN,xNS);
   end
else
   if((xN-n)/xNS > xNMN)
      p=0;return
   else
      p2=normpdf(n,xN,xNS);
   end
end
p = p1*p2;


