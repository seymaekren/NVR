
function [M,differenceMatrix] = NVR_RDC2PROB(TABLE,RDC,VECS,TEN,peakIndices,residueIndices);

%NVR_SHIFTX2PROB: This computes assignment probabilities based on RDCs, it is not meant
%             to be called by the user


%////////////////////////////////////////////////////////////////////////////////////////////
%//  NVR_RDC2PROB.m
%//
%//  Version:		0.1
%//
%//  Description:	 This computes assignment probabilities based on RDCs
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

%    NVR_RDC2PROB
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

%EPSILON = 1E-4;

%persistent firstRun

logicalM = TABLE;

[sizeX,sizeY]= size(logicalM);

M = ones(sizeX,sizeY) * 0.5;
for i = 1:sizeX
  for j = 1:sizeY
    M(i,j) = 1.0*logicalM(i,j);
  end
end


%M = TABLE;

%precompute rdcs for the tensor and vectors
backcomputed = zeros(1,size(VECS,1));
for(i=1:size(VECS,1))
   v = VECS(i,:);
   if (v(1) == -999)
     backcomputed(i) = -999;
   else
     backcomputed(i) = v*TEN*v';
   end
%   fprintf(1, 'i = %d vector: %f %f %f backcomputedRDC: %f\n', i, v(1),v(2),v(3),backcomputed(i));
end
%keyboard
% $$$ f1 = find(RDC>-999);
% $$$ if(range(RDC(f1))/range(backcomputed)<.9)
% $$$    backcomputed=backcomputed*(range(RDC(f1))/range(backcomputed));
% $$$ elseif(range(RDC(f1))/range(backcomputed)>1.1)
% $$$    backcomputed=backcomputed*(range(RDC(f1))/range(backcomputed));
% $$$ end

range(backcomputed)/10 ;
%keyboard
differenceMatrix = ones(size(M,1),size(M,2))*-1;%all unfilled
                                                %entries will be negative
for(i=1:size(M,1))
  
%  fprintf(1, 'i = %d\n',i);
  if(RDC(peakIndices(i))==-999)%means there was no RDC, give it a uniform probability over possible assignments
%      fprintf(1, 'skipping.\n');
      M(i,:) = 1.0*TABLE(i,:);
   else
   
   for(j=1:size(M,2))
      %compute the likelyhood of seeing the observed RDC given the vector 
         %and the tensor
         if(TABLE(i,j)>0)
	   %%%	   M(i,j) =
%normpdf(RDC(peakIndices(i))-backcomputed(residueIndices(j)),0,range(backcomputed)/8);
           if (backcomputed(residueIndices(j)) == -999)
	     M(i,j) = 0;
	     continue;
	   end


           M(i,j) = normpdf(RDC(peakIndices(i))-backcomputed(residueIndices(j)),0,1);
%          M(i,j) = normpdf(RDC(peakIndices(i))-backcomputed(residueIndices(j)),0,range(backcomputed)/8);
	  
%          if (isempty(firstRun))
%	    fprintf(1, 'setting the RDC std to range/8\n');
%	    keyboard
%	    firstRun = 0;
%	  end

	  
	  %	   if (M(i,j) == 0)
	  %	     M(i,j) = EPSILON;
	  %	     fprintf(1, 'warning. no longer setting 0 probs in RDC2PROB.\n');
	  %keyboard
	  %	   end
	   
	   if ((i == j) & (M(i,j) == 0))
	     fprintf(1, 'peak %d has 0 RDC prob.\n', i);
	     fprintf(1, 'RDC = %f backcomputed RDC = %f\n', ...
		     RDC(peakIndices(i)), ...
		     backcomputed(residueIndices(j)));
	   %  keyboard
% $$$ 	   elseif ((M(i,j) == 0))
% $$$ 	     fprintf(1, 'peak %d residue %d has 0 RDC prob.\n', i,j);
% $$$ 	     fprintf(1, 'RDC = %f backcomputed RDC = %f\n', ...
% $$$ 		     RDC(peakIndices(i)), ...
% $$$ 		     backcomputed(residueIndices(j)));
	   %  keyboard
	   end
	   
	   
	   
	   differenceMatrix(i,j) = 1/(1+exp(abs(RDC(peakIndices(i))-backcomputed(residueIndices(j)))));
% $$$ 	   if (i == j)
% $$$ 	     fprintf(1, 'i = %d, backcomputed RDC= %f\n', i, backcomputed(residueIndices(j)));
% $$$ 	   end
% $$$ 	   if (i == 24) & (j == 24)
% $$$ 	     fprintf(1, 'i=24, j=24\n');
% $$$ 	     fprintf(1, 'RDC_exp = %f RDC_bc = %f\n', RDC(peakIndices(i)), ...
% $$$ 		     backcomputed(residueIndices(j)));
% $$$ 	     fprintf(1, 'prob = %f\n', M(i,j));
% $$$ 	     keyboard
% $$$ 	   end
	   
	   
	   %	   if (M(i,j) == 0)
%	     fprintf(1, 'the backcomputed and exp. RDC are so far off.\n');
%	     fprintf(1, 'that the probability is 0.\n');
%	     fprintf(1, 'peakIndex = %d residueIndex = %d exp_rdc = %f bc_rdc = %f\n',peakIndices(i),residueIndices(j),RDC(peakIndices(i)),backcomputed(residueIndices(j)));
%	   end
	 end
   end
   end
   if(sum(M(i,:))>0)
      M(i,:) = M(i,:)/sum(M(i,:));%re-normalize   
   else
     fprintf(1, 'there is an RDC but the probabilities are 0.\n');
     fprintf(1, 'assigning to ASSIGNTABLE entries with equiprob.\n');
%     keyboard
     M(i,:) = 1.0*TABLE(i,:);
     if(sum(M(i,:))>0)
       M(i,:) = M(i,:)/sum(M(i,:));%re-normalize   
     else
       fprintf(1, 'both ASSIGNTABLE and RDC are giving 0 entries.\n');
       fprintf(1, 'assigning to equiprobability.\n');
       M(i,:) = 1/length(M(i,:));
     end
%     keyboard
   end
end

