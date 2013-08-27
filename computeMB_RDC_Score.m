
function [M,numRDC] = computeMB_RDC_Score(TABLE,RDC1,RDC2, VECS,TEN1, TEN2);

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


% $$$ logicalM = TABLE;
% $$$ 
% $$$ [sizeX,sizeY]= size(logicalM);
% $$$ 
% $$$ M = ones(sizeX,sizeY) * 0.5;
% $$$ for i = 1:sizeX
% $$$   for j = 1:sizeY
% $$$     M(i,j) = 1.0*logicalM(i,j);
% $$$   end
% $$$ end
% $$$ 

%M = TABLE;

%precompute rdcs for the tensor and vectors
backcomputed1 = zeros(1,size(VECS,1));
backcomputed2 = zeros(1,size(VECS,1));
for(i=1:size(VECS,1))
   v = VECS(i,:);
   backcomputed1(i) = v*TEN1*v';
   backcomputed2(i) = v*TEN2*v';
end

f1 = find(RDC1>-999);
if(range(RDC1(f1))/range(backcomputed1)<.9)
   backcomputed1=backcomputed1*(range(RDC1(f1))/range(backcomputed1));
elseif(range(RDC1(f1))/range(backcomputed1)>1.1)
   backcomputed1=backcomputed1*(range(RDC1(f1))/range(backcomputed1));
end

f2 = find(RDC2>-999);
if(range(RDC2(f2))/range(backcomputed2)<.9)
   backcomputed2=backcomputed2*(range(RDC2(f2))/range(backcomputed2));
elseif(range(RDC2(f2))/range(backcomputed2)>1.1)
   backcomputed2=backcomputed2*(range(RDC2(f2))/range(backcomputed2));
end

%Da1 = max(RDC1)/2; lowerThr1 = Da1 * 2.5/100; higherThr1 = Da1 * 0.1;
%Da2 = max(RDC2)/2; lowerThr2 = Da2 * 2.5/100; higherThr2 = Da2 *
%0.1;

Da1 = max(RDC1)/2; lowerThr1 = Da1 * 10/100; higherThr1 = Da1 * 40/100;
Da2 = max(RDC2)/2; lowerThr2 = Da2 * 10/100; higherThr2 = Da2 * 40/100;

[U1,Sigma1,V1] = svd(TEN1);
[U2,Sigma2,V2] = svd(TEN2);
Da1_prediction2 = max(diag(Sigma1))/2;
Da2_prediction2 = max(diag(Sigma2))/2;

%keyboard

M         = TABLE * 0;
%sigma_rdc = 0.21;
numRDC    = 0;

[numPeaks,numResidues] = size(M);
assert (numPeaks == length(RDC1));
assert (numResidues == size(VECS,1));


%fprintf(1, 'cancelled MB RDC Score part on RDC1.\n');
%keyboard
%if (0)

for(i=1:size(M,1))
   if (RDC1(i)==-999)%means there was no RDC, give it a uniform probability over possible assignments
     %M(i,:) = 1.0*TABLE(i,:);
   else
     numRDC = numRDC  + 1;
     for(j=1:size(M,2))
       %compute the likelyhood of seeing the observed RDC given the vector 
       %and the tensor
       if(TABLE(i,j)>0)

	 diff1 = abs(RDC1(i)-backcomputed1(j));
	 
	 if (diff1 < lowerThr1)
	   M(i,j) = M(i,j) + 1;
	 elseif (diff1 < higherThr1)
	   M(i,j) = (higherThr1-diff1) / (higherThr1-lowerThr1);
	 end
%	 keyboard

       end
     end
   end
%   if(sum(M(i,:))>0)
%      M(i,:) = M(i,:)/sum(M(i,:));%re-normalize   
%   end
end
%end

for(i=1:size(M,1))
   if (RDC2(i)==-999)%means there was no RDC, give it a uniform probability over possible assignments
     %M(i,:) = 1.0*TABLE(i,:);
   else
     numRDC = numRDC + 1;
     for(j=1:size(M,2))
       %compute the likelyhood of seeing the observed RDC given the vector 
       %and the tensor
       if(TABLE(i,j)>0)

	 diff2 = abs(RDC2(i)-backcomputed2(j));
	 
	 if (diff2 < lowerThr2)
	   M(i,j) = M(i,j) + 1;
	 elseif (diff2 < higherThr2)
	   M(i,j) = (higherThr2-diff2) / (higherThr2-lowerThr2);
	 end 
	 %	 keyboard
       end
     end
   end
   %   if(sum(M(i,:))>0)
   %      M(i,:) = M(i,:)/sum(M(i,:));%re-normalize   
   %   end
end
