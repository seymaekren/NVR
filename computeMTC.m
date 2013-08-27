function MTC = computeMTC(table, rdc1,rdc2, x, y, z, correctAssignments);

%computes the *median* tensor consistency.
%I think table contains the assignment table, rdc1 the experimental
%rdc values, x,y and z the internuclear bond vectors.

%computes two Saupe elements for a given assignment, by randomly 
%choosing two different mutually exclusive subsets of k elements.

% $$$ if (nargin < 6)
% $$$   numMaxTry = 1000;
% $$$ else
% $$$   numMaxTry;
% $$$ end


MTC1s      = computeMTCs(x,y,z,rdc1,table,correctAssignments);
MTC2s      = computeMTCs(x,y,z,rdc2,table,correctAssignments);
MTCs       = [MTC1s MTC2s];
%figure; 
%plot(MTCs,'*')
MTC        = median(MTCs);

function MTCs = computeMTCs(x,y,z,rdc1,table,correctAssignments)

maxNumIter = 100;
MTCs       = [];
[n,m]      = size(table);

S1 = 0;   S2 = 0;

S1 = updateTen(table,rdc1,[x,y,z]);

    
    if (numIter == maxNumIter)
      fprintf('error, cannot find a Saupe tensor matrix\n');
      keyboard;
    end
    
    numIter = numIter + 1;
    %      fprintf(1,'*');
    %    rp4 = randperm(n);
    %    selectedElements = rp4(1:k);
    
    selectedElements = rp3(k+1:2*k); 

    assert (2*k <= length(rp3));
    
    tableSubset2     = table(selectedElements, :);
    
    S2 = updateTen(tableSubset2,rdc1(selectedElements),[x,y,z]);
  end
  
  if (numIter == maxNumIter)
    fprintf('error, cannot find a Saupe tensor matrix\n');
    keyboard;
  end
  
  
% $$$   S2 = 0;
% $$$   numIter = 0;
% $$$   
% $$$   while ((S2 == 0) & (numIter < maxNumIter))
% $$$     numIter = numIter + 1;
% $$$     %      fprintf(1,'*');
% $$$     rp4 = randperm(n);
% $$$     selectedElements = rp4(1:k);
% $$$     
% $$$     tableSubset2 = table(selectedElements, :);
% $$$     
% $$$     S2 = updateTen(tableSubset2,rdc1(selectedElements),[x,y,z]);
% $$$   end
% $$$   
% $$$   if (numIter == maxNumIter)
% $$$     fprintf('error, cannot find a Saupe tensor matrix\n');
% $$$     keyboard;
% $$$   end
  
  
  assert (S2 ~= 0);
  MTCs = [MTCs real(NVR_COMP_TEN(S1, S2))];

end


%##################################################################################
%##  updateTen
function SAUPE = updateTen(TABLE,RDC,VECS)
SAUPE  = 0;
M = zeros(1,5);
vec = 0;
ct=1;
for(i=1:size(TABLE,1))
   x = find(TABLE(i,:));
   if(length(x)==1 & RDC(i)>-999)
      M(ct,1) = VECS(x,1).^2-VECS(x,3).^2;
      M(ct,2) = 2*(VECS(x,1)*VECS(x,2));
      M(ct,3) = 2*(VECS(x,1)*VECS(x,3));
      M(ct,4) = VECS(x,2).^2-VECS(x,3).^2;
      M(ct,5) = 2*(VECS(x,2)*VECS(x,3));
      vec(ct) =  RDC(i);
      ct=ct+1;
   end
end
if(size(M,1)<5)
   return
end
%invert the matrix
if(size(M,1)==5 & size(M,2)==5)
   M = inv(M);
else
   M = pinv(M);
end
%do the least squares fitting
s = M*vec';
SAUPE = zeros(3,3);
SAUPE(1,1) = s(1);
SAUPE(2,2) = s(4);
SAUPE(3,3) = -1*(s(1)+s(4));
SAUPE(1,2) = s(2);
SAUPE(2,1) = s(2);
SAUPE(1,3) = s(3);
SAUPE(3,1) = s(3);
SAUPE(2,3) = s(5);
SAUPE(3,2) = s(5);
