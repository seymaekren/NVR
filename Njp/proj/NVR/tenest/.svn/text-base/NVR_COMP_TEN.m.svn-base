 function percentile = NVR_COMP_TEN(Saupe1,Saupe2)
 
 %NVR_COMP_TEN: This returns the percentile reporting the similarity of two tensors. Higher values indicate
 %					 higher degrees of similarity. 
%      Input:  Saupe1- a 3x3 saupe matrix
%					Saupe2- a 3x3 saupe matrix
%		 Output: percentile - The  the percentile reporting the similarity of two tensors.  


%////////////////////////////////////////////////////////////////////////////////////////////
%//  NVR_COMP_TEN.m
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

%    NVR_COMP_TEN

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

 
 eigval1  = zeros(3,3);
eigvect1 = zeros(3,3);
[eigvect1,eigval1] = eig(Saupe1);
[dg1,ind]=sort(diag(eigval1));


Vx1 = eigvect1(1:3,ind(2));
Vx1 = Vx1/norm(Vx1);

Vy1 = eigvect1(1:3,ind(1));
Vy1 = Vy1/norm(Vy1);

Vz1 = eigvect1(1:3,ind(3));
Vz1 = Vz1/norm(Vz1);

   
eigval2  = zeros(3,3);
eigvect2 = zeros(3,3);
[eigvect2,eigval2] = eig(Saupe2);
[dg1,ind]=sort(diag(eigval2));


Vx2 = eigvect2(1:3,ind(2));
Vx2 = Vx2/norm(Vx2);

Vy2 = eigvect2(1:3,ind(1));
Vy2 = Vy2/norm(Vy2);

Vz2 = eigvect2(1:3,ind(3));
Vz2 = Vz2/norm(Vz2);



AngleX = acos(dot(Vx1,Vx2));
AngleY = acos(dot(Vy1,Vy2));
AngleZ = acos(dot(Vz1,Vz2));

AngleX  = min(AngleX , pi-AngleX);
AngleY  = min(AngleY , pi-AngleY);
AngleZ  = min(AngleZ , pi-AngleZ);


Angles = sort([AngleX AngleY AngleZ]);

SXX=AngleX*(180/pi);
SYY=AngleY*(180/pi);
SZZ=AngleZ*(180/pi);

alpha = Angles(3);
beta  = Angles(2);

prob1 = (2*pi)*(1-cos(alpha))/(4*pi);

prob2 = 2*beta/(2*pi);

degen = 4;

prob = degen*prob1*prob2;

percentile = 1-prob;

