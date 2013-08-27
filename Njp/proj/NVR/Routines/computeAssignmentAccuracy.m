function [assignmentAccuracy,assignments] = computeAssignmentAccuracy(peakIDs, RESNUMS, MASTER)

a = load('answerkey.m');
assignments=0;
correct=0;
incorr=0;

for(i=1:size(MASTER,1))

  pk               =peakIDs(i);  %get the peak id
  assignments(i,1) =pk;
   
  x = find(MASTER(i,:));
   
  if (isempty(x))
    rn = -1;      %peak i is unassigned.
  else
    rn=RESNUMS(x);%get the residue it was assigned to.
  end

  assignments(i,2)=rn;
  
  foo = find(a(:,1)==pk);
  if(rn==a(foo,2))
    correct=correct+1;
    incorr(i)=0;
  else
    incorr(i)=1;
  end
end

assignmentAccuracy = correct/size(MASTER,1)*100;

printAssignmentAccuracy = 1;
if (printAssignmentAccuracy)

  fprintf('\n');
  fprintf('\n');
  fprintf('Assignment Accuracy = %f percent \n',correct/size(MASTER,1)*100);
  fprintf('\n');
  fprintf('\n');
  fprintf('Assignments\n');
  fprintf('Peak ID -> Residue Number\n');
  for(i=1:size(assignments,1))
    if(incorr(i)==0)
      fprintf('%d	%d\n',assignments(i,1),assignments(i,2));
    else
      fprintf('*%d	%d  (correct=%d->%d)\n',assignments(i,1),assignments(i,2),a(i,1),a(i,2));
    end
  end
  fprintf('\n');
  fprintf('\n');
  %%keyboard
  if(size(MASTER,1)~=size(MASTER,2))
    %find which residues aren't assigned
    ct=1;
    missing=[];
    
    for(i=1:length(RESNUMS))
      rn = RESNUMS(i);
      if(length(find(assignments(:,2)==rn))==0)
	missing(ct)=rn;
	ct=ct+1;
      end
    end
    
    fprintf('The following residues are missing peaks\n');
    for(i=1:length(missing))
      fprintf('%d\n',missing(i));
    end
  end
  
end
