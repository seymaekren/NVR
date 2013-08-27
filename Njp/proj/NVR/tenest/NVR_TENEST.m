function [S]=NVR_TENEST(rdcs,vectors)

%NVR_TENEST: This program estimates the alignment tensor for a given medium using 
%      Input:  rdcs- a 1xN matrix containing N rdcs from a single medium, in an arbitrary order.
%					vectors- a Mx3 matrix containing M 1 x 3 normalized vectors corresponding to
%					the bond vectors (taken from a PDB model). The order of the vectors is arbitrarty, 
%					but the vectors should, of course, correspond to the actual spin systems associated
%					with the RDC data. 
%		 Output: S- The 3x3 estimated tensor.  


%////////////////////////////////////////////////////////////////////////////////////////////
%//  NVR_TENEST.m
%//
%//  Version:		0.1
%//
%//  Description:	 estimates the alignment tensor for a given medium using bond vectors from a model
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

%    NVR_TENEST

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


rdcs = rdcs(find(rdcs>-999));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%First, estimate Da and Dr using the powder pattern method
[maxdz,maxdy] = pp(rdcs, size(vectors,1));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Next, estimate the eigenvectors

S = KLSEARCH(maxdz,maxdy, rdcs, vectors);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function BSAUPE = KLSEARCH(DZ,DY, rdcs, vectors)

%construct eigenvalue matrix
EIGS = zeros(3,3);
EIGS(1,1)=DY;
EIGS(2,2)=-1*(DZ+DY);
EIGS(3,3)=DZ;

BESTKL=1000000000;



%this is a tunable parameter, basically, it sets the number of bins
%in the histogram which is used to generate the probability distributions
HISTSIZE=15;

%construct reference histogram
ohist = hist(rdcs,HISTSIZE);
ohist = ohist/sum(ohist);


%this is a tunable parameter, basically, it determines how finely you sample the sphere
%now, generate rotations 
numsamps=36;

%this is a tunable parameter, basically, it determines how finely you sample rotations
%its the number of parts you divide 360 degrees into
thetres = 36;


[x,y,z] = NVR_sph2(numsamps);



%try all rotations, select the best one
for(i=1:size(x,2))
   
   %construct the coordinate system
   xaxis= normalize([x(i) y(i) z(i)]);
   yaxis = normalize(cross([0 0 1],xaxis)); 
   zaxis = normalize(cross(xaxis,yaxis));
   
   %we now have our coordinate system, try all rotations, 10 degree increments
   for(theta=0:(2*pi)/thetres:2*pi)
      skip=0;
      
      %construct rotation matrix
      x1 = xaxis;
      y1 = cos(theta)*yaxis + sin(theta)*zaxis;
      z1 = -1*sin(theta)*yaxis + cos(theta)*zaxis;
      ROT = [y1;x1;z1]';
      
      %now, for each vector in the model , rotate it, and compute expected rdc
      predrdcs = zeros(1,size(vectors,1));
      SAUPE= ROT*EIGS*ROT';
      for(j=1:size(vectors,1))
         v = vectors(j,:);
         predrdcs(j) = v * SAUPE *v';
      end      
      
      %construct putative histogram
      phist = hist(predrdcs,HISTSIZE);
      phist = phist/sum(phist);
      
      %compute the relative entropy
      kldist = KL(ohist,phist);
      
      if(kldist < BESTKL)
         BESTKL = kldist ;  
         BSAUPE=SAUPE;
      end
   end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [maxdz,maxdy] = pp(rdcs,proteinsize);
maxscore=-1000000000;

maxdz =-999;
maxdy =-999;
%these are adjustable parameters. Basically, they define the 
%window about the max and min observed RDC value overwhich we 
%search for the 'true' max and min value. These, in turn, define
%the rhombicity of the tensor
maxlo = 2;
maxhi = 2;
minlo = 2;
minhi = 2;


for(dz = max(rdcs)-maxlo:.5:max(rdcs)+maxhi)
   for(dy = min(rdcs)-minlo :.5:min(rdcs)+minhi)
      
      dx = -1*(dz+dy);
      da = (1/2)*dz;
      dr = (1/3)*(dx-dy);
      R = dr/da;
      
      samp=size(rdcs,1);
      probs = 0;
      ct =1;
      
      
      if(dz>dx & dx >=dy)
         for(i=1:size(rdcs,1))
            
            %type 1 soln
            if(rdcs(i) > dx & rdcs(i)<=dz & rdcs(i)>=dy)
               
               xisq = (2*da -rdcs(i))/(3*da-(3/2)*dr);
               etasq = (2*da -rdcs(i))/(3*da+(3/2)*dr);
               lamsq = (rdcs(i)+da+(3/2)*dr)/(3*dr);
               zetsq = (rdcs(i)+da+(3/2)*dr)/(3*da+(3/2)*dr); 
               
               xi = sqrt(xisq);
               eta = sqrt(etasq);
               lam = sqrt(lamsq);
               zet = sqrt(zetsq);
               
               A = 2/(3*pi*da*zet);
               B = sqrt(1/(4-R^2));
               k = xi/lam;
               m = k^2;
               if(m>=0 & m<=1 )
                  C = ellipke(m);
                  probs(ct) =A*B*C;  %compute likelihood
                  if(isinf(probs(ct))==1)
                  else
                     ct = ct+1;
                  end
               end
               %type 2 soln   
            elseif(rdcs(i)<=dz & rdcs(i)>=dy)
               xisq = (2*da -rdcs(i))/(3*da-(3/2)*dr);
               etasq = (2*da -rdcs(i))/(3*da+(3/2)*dr);
               lamsq = (rdcs(i)+da+(3/2)*dr)/(3*dr);
               zetsq = (rdcs(i)+da+(3/2)*dr)/(3*da+(3/2)*dr);
               
               xi = sqrt(xisq);
               eta = sqrt(etasq);
               lam = sqrt(lamsq);
               zet = sqrt(zetsq);
               
               A = 1/(3*pi*da*eta);
               B = sqrt(2/(2*R+R^2));
               k = lam/xi;
               m = k^2;
               if(m>=0 & m<=1)
                  C = ellipke(m);
                  probs(ct) =A*B*C; %compute likelihood
                  if(isinf(probs(ct))==1)
                  else
                     ct = ct+1;
                  end
               end
            end
         end
         
         %normalize the length of the vector of probabilities 
         %suth that its the same length as the size of the protein
         s = 0;
         if(length(probs)<proteinsize)
            minsc = min(probs)/10;  %give the 'missing ones' some arbitrary small value
            for(ct=length(probs)+1:proteinsize)
               probs(ct)=minsc;
            end
         end
         
         %get total log-likelihood
         for(i=1:length(probs))
            s = s+log(probs(i));
         end
         
         %keep track of the best ones
         if(s>maxscore)
            maxscore=s;
            maxdz = dz;
            maxdy = dy;
         end
         
      end
   end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = KL(d1,d2)

%re-normalize
d1=d1+.00000001;
d2=d2+.00000001;

d1 = d1/sum(d1);
d2 = d2/sum(d2);

s1 = 0; 
s2 = 0;

for(i=1:length(d1))
   s1 = s1+d1(i)*(log(d1(i)/d2(i)));
   s2 = s2+d2(i)*(log(d2(i)/d1(i)));
end

s = (abs(s1)+abs(s2));

