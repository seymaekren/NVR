function M = HD_TOCSY2PROB(peakIDs,H,N,TYPES,SSTRUCT,NOES,ALLDISTS,NTH,ROWIN,COLIN);

%NVR_TOCSY2PROB: This computes assignment probabilities based on BMRB side-chain statistics, it is not meant
%             to be called by the user


%////////////////////////////////////////////////////////////////////////////////////////////
%//  NVR_TOCSY2PROB.m
%//
%//  Version:		0.1
%//
%//  Description:	 This computes assignment probabilities based on BMRB statistics
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

%    NVR_TOCSY2PROB
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


[TH1,TH2,TN,TRN]=textread('TOCSY.m','%f %f %f %f');

M = zeros(length(H),length(TYPES));


for(i=1:size(M,1))
   
   %get the number of tocsy peaks
   rn=peakIDs(i);
   NUMTOCSYPEAKS=length(find(TRN==rn));
   TOCSYPEAKS=TH2(find(TRN==rn));
   
   
   for(j=1:length(TYPES))
      if(strcmp(TYPES(j),'ALA')==1)
         if(NUMTOCSYPEAKS>3 )
            M(i,j)=0; 
         else
            M(i,j)=getTOCSYProb(TOCSYPEAKS,'A',SSTRUCT(j),length(TYPES));
         end   
      elseif(strcmp(TYPES(j),'CYS')==1)
         if(NUMTOCSYPEAKS>5)
            M(i,j)=0; 
         else
            M(i,j)=getTOCSYProb(TOCSYPEAKS,'C',SSTRUCT(j),length(TYPES));
         end   
      elseif(strcmp(TYPES(j),'ASP')==1)
         if(NUMTOCSYPEAKS>4)
            M(i,j)=0; 
         else
            M(i,j)=getTOCSYProb(TOCSYPEAKS,'D',SSTRUCT(j),length(TYPES));
         end
      elseif(strcmp(TYPES(j),'GLU')==1)
         if(NUMTOCSYPEAKS>6)
            M(i,j)=0; 
         else
            M(i,j)=getTOCSYProb(TOCSYPEAKS,'E',SSTRUCT(j),length(TYPES));
         end   
      elseif(strcmp(TYPES(j),'PHE')==1)
         
         if(NUMTOCSYPEAKS>9)
            M(i,j)=0; 
         else
            M(i,j)=getTOCSYProb(TOCSYPEAKS,'F',SSTRUCT(j),length(TYPES));
         end   
      elseif(strcmp(TYPES(j),'GLY')==1)
         if(NUMTOCSYPEAKS>3)
            M(i,j)=0; 
         else
            M(i,j)=getTOCSYProb(TOCSYPEAKS,'G',SSTRUCT(j),length(TYPES));
         end
      elseif(strcmp(TYPES(j),'HIS')==1)
         if(NUMTOCSYPEAKS>8)
            M(i,j)=0; 
         else
            M(i,j)=getTOCSYProb(TOCSYPEAKS,'H',SSTRUCT(j),length(TYPES));
         end
      elseif(strcmp(TYPES(j),'ILE')==1)
         if(NUMTOCSYPEAKS>7)
            M(i,j)=0; 
         else
            M(i,j)=getTOCSYProb(TOCSYPEAKS,'I',SSTRUCT(j),length(TYPES));
         end
      elseif(strcmp(TYPES(j),'LYS')==1)
         if(NUMTOCSYPEAKS>11)
            M(i,j)=0; 
         else
            M(i,j)=getTOCSYProb(TOCSYPEAKS,'K',SSTRUCT(j),length(TYPES));
         end      
      elseif(strcmp(TYPES(j),'LEU')==1)
         if(NUMTOCSYPEAKS>7)
            M(i,j)=0; 
         else
            M(i,j)=getTOCSYProb(TOCSYPEAKS,'L',SSTRUCT(j),length(TYPES));
         end
      elseif(strcmp(TYPES(j),'MET')==1)
         if(NUMTOCSYPEAKS>7)
            M(i,j)=0; 
         else
            M(i,j)=getTOCSYProb(TOCSYPEAKS,'M',SSTRUCT(j),length(TYPES));
         end
      elseif(strcmp(TYPES(j),'ASN')==1)
         if(NUMTOCSYPEAKS>6)
            M(i,j)=0; 
         else
            M(i,j)=getTOCSYProb(TOCSYPEAKS,'N',SSTRUCT(j),length(TYPES));
         end
      elseif(strcmp(TYPES(j),'GLN')==1)
         if(NUMTOCSYPEAKS>8)
            M(i,j)=0; 
         else
            M(i,j)=getTOCSYProb(TOCSYPEAKS,'Q',SSTRUCT(j),length(TYPES));
         end
      elseif(strcmp(TYPES(j),'ARG')==1)
         M(i,j)=getTOCSYProb(TOCSYPEAKS,'R',SSTRUCT(j),length(TYPES));
      elseif(strcmp(TYPES(j),'SER')==1)
         if(NUMTOCSYPEAKS>5)
            M(i,j)=0; 
         else
            M(i,j)=getTOCSYProb(TOCSYPEAKS,'S',SSTRUCT(j),length(TYPES));
         end
      elseif(strcmp(TYPES(j),'THR')==1)
         if(NUMTOCSYPEAKS>5)
            M(i,j)=0; 
         else
            M(i,j)=getTOCSYProb(TOCSYPEAKS,'T',SSTRUCT(j),length(TYPES));
         end
      elseif(strcmp(TYPES(j),'VAL')==1)
         if(NUMTOCSYPEAKS>5)
            M(i,j)=0; 
         else
            M(i,j)=getTOCSYProb(TOCSYPEAKS,'V',SSTRUCT(j),length(TYPES));
         end
      elseif(strcmp(TYPES(j),'TRP')==1)
         if(NUMTOCSYPEAKS>10)
            M(i,j)=0; 
         else
            M(i,j)=getTOCSYProb(TOCSYPEAKS,'W',SSTRUCT(j),length(TYPES));
         end
      elseif(strcmp(TYPES(j),'TYR')==1)
         if(NUMTOCSYPEAKS>9)
            M(i,j)=0; 
         else
            M(i,j)=getTOCSYProb(TOCSYPEAKS,'Y',SSTRUCT(j),length(TYPES));
         end
      else
         PROBLEM = TYPES(j)   
      end
      
   end
   
   M(i,:) = M(i,:)/sum(M(i,:));%re-normalize
end

M = thresh(M,min(nonzeros(M))*10e+40);

TABLE = and(M,M);
nlast = sum(sum(TABLE));
for(i=1:100)
   NP = NVR_NOE2PROB(TABLE(1:size(M,1),:),NOES,ALLDISTS,NTH,ROWIN,COLIN);
   %note that noe pruning is internally called here.
   TABLE(1:size(M,1),:)=and(TABLE(1:size(M,1),:),NP);
   if(sum(sum(TABLE)) == nlast)
      break;
   end
   nlast = sum(sum(TABLE));
   M = M.*TABLE;
   for(i=1:size(M,1))
      M(i,:)=M(i,:)/sum(M(i,:));
   end
end

%renornmalize
for(i=1:size(M,1))
   M(i,:)=M(i,:)/sum(M(i,:));
end



function tp=getTOCSYProb(TOCSYPEAKS,AATYPE,SSTYPE,len)
nm=sprintf('/home/home4/apaydin/Mist/NVR/HDB/%s%s.mat',char(AATYPE),char(SSTYPE));

SHIFTS=load(nm);SHIFTS=SHIFTS.shifts;
MB=zeros(length(SHIFTS(1,:)))+999;
for(i=1:length(TOCSYPEAKS))
   for(j=1:size(SHIFTS,2))
      MB(i,j)=abs(TOCSYPEAKS(i)-SHIFTS(1,j));
   end
end

h=hungarian(MB');
fprintf(1, 'size(MB,1) = %d, size(MB,2) = %d\n',size(MB,1),size(MB,2));
keyboard
tp=1;
for(i=1:length(TOCSYPEAKS))
   pos=h(i);
   if(MB(i,pos)<100)
      if(SHIFTS(2,pos)==0)
         p=0;
      else
         p=MB(i,pos)/SHIFTS(2,pos);
      end
      tp=tp*(2*(1-tcdf(p,len)));%i should probably have just used pdf here
   end
end


