function M = NVR_HD2PROB(TABLE,HDEXCHANGE,HDBOND);

%NVR_HD2PROB: This computes assignment probabilities based on amide exchange rates, it is not meant
%             to be called by the user


%////////////////////////////////////////////////////////////////////////////////////////////
%//  NVR_HD2PROB.m
%//
%//  Version:		0.1
%//
%//  Description:	 This computes assignment probabilities based on amide exchange rates
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

%    NVR_HD2PROB
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

M=and(TABLE,TABLE);
for(i=1:size(TABLE,1))
   if(HDEXCHANGE(i)==1)
      for(j=1:length(HDBOND))
         if(strcmp(HDBOND(j),'N')==1)
            M(i,j)=0; 
         end
      end
   end
   if(sum(M(i,:))==0)
      M(i,:)=1;
   end
   M(i,:) = M(i,:)/sum(M(i,:));%re-normalize
end
