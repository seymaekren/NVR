function [x,y,z] = NVR_sph2(n)
%NVR_sph2: This program that (approximately) uniformly samples the unit sphere  
%      Input:  n- the number of points to sample on the sphere
%		 Output: [x,y,z]- The x,y, and z coordinates of the samples


%////////////////////////////////////////////////////////////////////////////////////////////
%//  NVR_sph2.m
%//
%//  Version:		0.1
%//
%//  Description:	 This program that (approximately) uniformly samples the unit sphere  
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

%    NVR_sph2
%    Copyright (C) 2002  Christopher James Langmead

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
% k is the number of points around the equator



k=round(sqrt(pi*n));

% Compute angle deltas
dt = 2*pi/k;

x=zeros(1,1);
y=zeros(1,1);
z=zeros(1,1);

c=1;
phi=0;
theta=dt;
for i=1:round(k/2-1)
   phi=0;
   pts = round(k*sin(theta));
   dp = 2*pi/pts;
   for j=1:pts
      x(c)=sin(theta)*cos(phi);
      y(c)=sin(theta)*sin(phi);
      z(c)=cos(theta);
      phi=phi+dp;
      c=c+1;      
   end
   theta=theta+dt;
end

