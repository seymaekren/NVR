function MASTER = checkMASTER(MASTER,S1,RDC1,VECTORS,S2,RDC2,TP_A,CP_A,SXCP_A,SSCP_A,HDE_A,RP1_A,RP2_A);

%checks whether the MASTER matrix contains more than one assigned
%peak or residue. If so, it recomputes the assignment only for
%those peaks and residues which have multiple assignments.

if(max(max(sum(MASTER)), max(sum(MASTER')))>1)
  for(i=1:size(MASTER,1))
       if(length(find(MASTER(i,:)))>1)
          MASTER(i,:)=0;
       end
  end
  
  for(i=1:size(MASTER,2))
    if(length(find(MASTER(:,i)))>1)
      MASTER(:,i)=0;
    end
  end
  
  S1 = updateTen(MASTER,RDC1,VECTORS);
  S2 = updateTen(MASTER,RDC2,VECTORS);
  
  cand = MASTER(1,:)*0;
  %now, find unassigned residues
  for(i=1:size(MASTER,2))
    if(length(find(MASTER(:,i)))==0)
      cand(i)=1;
    end
  end
  for(i=1:size(MASTER,1))
    if(length(find(MASTER(i,:)))==0)
      MASTER(i,:)=cand;
    end
  end
  
  X = TP_A.*CP_A.*SXCP_A.*SSCP_A.*HDE_A.*RP1_A.*RP2_A.*MASTER;
  h = hungarian(makesquare(X)*-1);
  MASTER = MASTER *0;
  for(i=1:size(MASTER,1))
    MASTER(i,h(i))=1;
  end
end
