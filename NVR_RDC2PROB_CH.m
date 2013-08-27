
function M = NVR_RDC2PROB_CH(TABLE,RDC,VECS,TEN,peakIndices,residueIndices);

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
   backcomputed(i) = v*TEN*v';
end

f1 = find(RDC>-999);
if(range(RDC(f1))/range(backcomputed)<.9)
   backcomputed=backcomputed*(range(RDC(f1))/range(backcomputed));
elseif(range(RDC(f1))/range(backcomputed)>1.1)
   backcomputed=backcomputed*(range(RDC(f1))/range(backcomputed));
end

range(backcomputed)/10 ;
%keyboard
for(i=1:size(M,1))
   if(RDC(peakIndices(i))==-999)%means there was no RDC, give it a uniform probability over possible assignments
      M(i,:) = 1.0*TABLE(i,:);
   else
   
   for(j=1:size(M,2))
      %compute the likelyhood of seeing the observed RDC given the vector 
         %and the tensor
         if(TABLE(i,j)>0)
	   if (i == 56) & ( j == 56)
% $$$ 	     RDC(peakIndices(i))
% $$$ 	     residueIndices(j)
%	     keyboard;
	   end
%%%	   M(i,j) =
%normpdf(RDC(peakIndices(i))-backcomputed(residueIndices(j)),0,range(backcomputed)/8);
	   M(i,j) = normpdf(RDC(peakIndices(i))-backcomputed(residueIndices(j)),0,4);
         end
      end
   end
   if(sum(M(i,:))>0)
      M(i,:) = M(i,:)/sum(M(i,:));%re-normalize   
   end
end

