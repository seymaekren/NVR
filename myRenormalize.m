function [BPG_TO_UPDATE, impossibleToAssign] = myRenormalize(BPG_TO_UPDATE)

impossibleToAssign = 0;

for (i=1:size(BPG_TO_UPDATE,1))
  if (sum(BPG_TO_UPDATE(i,:))==0)
%    fprintf(1, 'row of bpg is empty.\n');
%    impossibleToAssign = 1;
%    keyboard
%    return;
  else
    BPG_TO_UPDATE(i,:)=BPG_TO_UPDATE(i,:) / sum(BPG_TO_UPDATE(i,:));
  end
end