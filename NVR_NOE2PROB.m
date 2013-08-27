function [possibleAssignmentsBPG, impossibleToAssign] = NVR_NOE2PROB(possibleAssignmentsBPG,NOES,ALLDISTS,NTH, ROWIN,COLIN);


%NVR_NOE2PROB: This computes assignment probabilities based on NOES


%////////////////////////////////////////////////////////////////////////////////////////////
%//  NVR_NOE2PROB.m
%//
%//  Version:		0.1
%//
%//  Description:	 This computes assignment probabilities based on noes
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

%    NVR_NOE2PROB
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


impossibleToAssign = 0;
possibleAssignmentsBPG=and(possibleAssignmentsBPG,possibleAssignmentsBPG);

for relPeakIndex=1:size(possibleAssignmentsBPG,1)

  absolutePeakIndicesOfPeaksHavingAnNOEWithPeak_relPeakIndex = find(NOES(ROWIN(relPeakIndex),:));%set of peaks having an noe with this one
  
  if (length(absolutePeakIndicesOfPeaksHavingAnNOEWithPeak_relPeakIndex) ...
      == 0 )      %if (peak # relPeakIndex has NOE constraints)
		  %NOES(ROWIN(relPeakIndex),:) refers
		  %to the NOES peak #relPeakIndex has.
		  continue;
  end
    
  possibleAbsoluteResidueIndicesForPeak_RelPeakIndex = COLIN(find(possibleAssignmentsBPG(relPeakIndex,:))); %find out who it is could be
          
  if (isempty(possibleAbsoluteResidueIndicesForPeak_RelPeakIndex))
%    impossibleToAssign = 1;
%    fprintf(1, 'found impossibleToAssign in NVR_NOE2PROB\n');
%    keyboard SB. This code could be different for pure a*
%    return;
%    fprintf(1, 'found empty ASSIGNTABLE row (for peak %d). continuing.\n',relPeakIndex);
%    keyboard
    continue;
  end

  for(j=1:length(absolutePeakIndicesOfPeaksHavingAnNOEWithPeak_relPeakIndex)) 
    %go through all NOE peaks with peak relPeakIndex
    
    relPeakIndexHavingNOEWithPeak_relPeakIndex=find(ROWIN==absolutePeakIndicesOfPeaksHavingAnNOEWithPeak_relPeakIndex(j));
    
    %previously it may have been possible not to find the peak having
    %an NOE with peak relPeakIndex, since that peak could have been
    %assigned. In that case the old HD did not make use of that
    %information. This allowed not to make incorrect pruning in case
    %previous assignment was incorrect, but also did not take advantage
    %of extra data, in case the previous assignment was correct.
    
    if(length(relPeakIndexHavingNOEWithPeak_relPeakIndex) == 0)
      %if peak with relative index
      %relPeakIndexHavingNOEWithPeak_relPeakIndex is assigned, then it
      %may be possible not to find its index in ROWIN.
      continue;
    end
    
%    if (relPeakIndex == 7) & (relPeakIndexHavingNOEWithPeak_relPeakIndex ...
%			      == 8)
%      fprintf(1, 'this case must be pruned.\n');
      %keyboard
%    end
    
    
    %            szin=length(find(M(relPeakIndexHavingNOEWithPeak_relPeakIndex,:)));
      
    relResidueIndices = find(possibleAssignmentsBPG(relPeakIndexHavingNOEWithPeak_relPeakIndex,:));
      
    %	  for (relResidueIndex=1:size(possibleAssignmentsBPG,2))
    for k = 1:length(relResidueIndices)
      relResidueIndex = relResidueIndices(k);
      
      assert (possibleAssignmentsBPG(relPeakIndexHavingNOEWithPeak_relPeakIndex,relResidueIndex)~= 0);
      
      absResidueIndex = COLIN(relResidueIndex);
      
      distOfRes_AbsResIdxToPotOfPk_RelPkIdx = ALLDISTS(possibleAbsoluteResidueIndicesForPeak_RelPeakIndex,absResidueIndex);
      
      if(min(distOfRes_AbsResIdxToPotOfPk_RelPkIdx)>NTH)
	%if residue # absoluteResidueIndex is too far away from potential
	%assignments to peak # relPeakIndex...
	%	      if (prunedPossibleAssignmentsBPG(relPeakIndexHavingNOEWithPeak_relPeakIndex,relResidueIndex)~= 0)

% $$$        	fprintf(1, 'peak %d has noe with peak %d\n', ROWIN(relPeakIndex),absolutePeakIndicesOfPeaksHavingAnNOEWithPeak_relPeakIndex(j));
% $$$        	fprintf(1, 'and the min. distance of the residues to which the peak can be assigned(');
% $$$        	for (myIndex = 1:length(possibleAbsoluteResidueIndicesForPeak_RelPeakIndex))
% $$$        	  fprintf(' %d, ',possibleAbsoluteResidueIndicesForPeak_RelPeakIndex(myIndex));
% $$$        	end
% $$$        	fprintf(1, ')to the residue (%d) is %f\n',absResidueIndex,min(distOfRes_AbsResIdxToPotOfPk_RelPkIdx));
% $$$        	fprintf(1, 'so pruning the assignment of %d to %d\n', absolutePeakIndicesOfPeaksHavingAnNOEWithPeak_relPeakIndex(j), absResidueIndex);
% $$$ 
% $$$ 	if (absolutePeakIndicesOfPeaksHavingAnNOEWithPeak_relPeakIndex(j) ...
% $$$ 	    == 13)
% $$$ 	  keyboard
% $$$ 	end
	
	possibleAssignmentsBPG(relPeakIndexHavingNOEWithPeak_relPeakIndex,relResidueIndex)=0;

	
	
	if (k == length(relResidueIndices)) %not checking before
                                            %trying the last residue
	  possibleAssignments = find(possibleAssignmentsBPG(relPeakIndexHavingNOEWithPeak_relPeakIndex,:));
	  if (isempty(possibleAssignments))
%	    fprintf(1,'pruned all assignments. impossibleToAssign is1\n');
	    impossibleToAssign = 1;
%	    keyboard
	    return;
	  end
	end
%	keyboard
      end %if
	  %            szout=length(find(M(relPeakIndexHavingNOEWithPeak_relPeakIndex,:)));
    end %for k
  end %for j
end%for relPeakIndex






