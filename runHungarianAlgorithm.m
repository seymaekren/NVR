%assignmentType is 1 for making assignments. it is 2 to retain the
%top top candidates for each each row.
function assignmentsBipartiteGraph = runHungarianAlgorithm(bipartiteGraphToVoteOn)
% $$$ if(size(ASSIGNTABLE,1)-size(ASSIGNTABLE,2)~=0)
% $$$ 
% $$$   ASSIGNTABLE = makeASSIGNTABLE_SquareAndNormalized(ASSIGNTABLE);
% $$$ 
% $$$ end 
% $$$ 
% $$$ ASSIGNTABLE=and(ASSIGNTABLE,ASSIGNTABLE);

if(size(bipartiteGraphToVoteOn,1)-size(bipartiteGraphToVoteOn,2)~=0)

  assert (size(bipartiteGraphToVoteOn,1) < size(bipartiteGraphToVoteOn,2)); ...
      %more residues than peaks
  
  squareBipartiteGraphToVoteOn =   makeBipartiteGraphToVoteOnSquare(bipartiteGraphToVoteOn);
  [squareBipartiteGraphToVoteOn, impossibleToAssign] =   myRenormalize(squareBipartiteGraphToVoteOn);

  assert (impossibleToAssign == 0);
  
else
  squareBipartiteGraphToVoteOn =   bipartiteGraphToVoteOn;
end

assignmentsBipartiteGraph      =   squareBipartiteGraphToVoteOn*0;


[assignmentVector,T]  = hungarian(squareBipartiteGraphToVoteOn'*-1);
[assignmentVector2,T2] = hungarian(squareBipartiteGraphToVoteOn*-1);

EPSILON = 1E-6;
assert (abs(T - T2) < EPSILON);
%for peakIndex=1:size(bipartiteGraphToVoteOn,1)
%  assignedResidueIndex = assignmentVector(peakIndex);
%  assert (assignmentVector2(assignedResidueIndex) == peakIndex);
%end

cost = 0;
for peakIndex=1:size(squareBipartiteGraphToVoteOn,1)
  cost = cost - squareBipartiteGraphToVoteOn(peakIndex,assignmentVector(peakIndex));
end
assert (abs(cost - T) < EPSILON);

cost2 = 0;
for residueIndex=1:size(squareBipartiteGraphToVoteOn,2)
  cost2 = cost2 - squareBipartiteGraphToVoteOn(assignmentVector2(residueIndex),residueIndex);
end
assert (abs(cost2 - T2) < EPSILON);




%fprintf(1, 'assert in hungarianrun passed.\n');
%keyboard
% $$$ elseif(assignmentType==2)
% $$$ 
% $$$   fprintf(1, 'currently, i dont think code is ever executed\n');
% $$$   keyboard
% $$$   
% $$$   sz=top;
% $$$    for(i=1:size(HM,1))
% $$$       [x ind]= sort(HM(i,:));
% $$$       in2=find(x);
% $$$       ind=ind(in2);x=x(in2);
% $$$       if(range(x)>0)
% $$$          if(length(x)<sz+1)
% $$$             assignments(i,ind)=1;
% $$$          else
% $$$             assignments(i,ind(length(ind)-sz:length(ind)))=1;
% $$$          end
% $$$       else
% $$$          assignments(i,ind)=1;
% $$$       end
% $$$       assignmentVector=[];   
% $$$    end
% $$$ else
% $$$   fprintf(1, 'I dont think this part of the code is reached\n');
% $$$   keyboard;
% $$$ %  assignmentVector=simp(HM);
% $$$ end
for(i=1:length(assignmentVector))
   assignmentsBipartiteGraph(i,assignmentVector(i)) = 1;
end

