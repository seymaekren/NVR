function [squareBipartiteGraphToVoteOn] = makeBipartiteGraphToVoteOnSquare (bipartiteGraphToVoteOn)

F=ones(max(size(bipartiteGraphToVoteOn)));
F(1:size(bipartiteGraphToVoteOn,1),1:size(bipartiteGraphToVoteOn,2))=bipartiteGraphToVoteOn;
squareBipartiteGraphToVoteOn=F;   


% $$$ 
% $$$ for(i=1:size(squareBipartiteGraphToVoteOn,1))
% $$$   if(sum(squareBipartiteGraphToVoteOn(i,:))==0)
% $$$     %i don't think the following line is reached but is ok to keep it.
% $$$     squareBipartiteGraphToVoteOn(i,:)=1;
% $$$   end
% $$$   squareBipartiteGraphToVoteOn(i,:)=squareBipartiteGraphToVoteOn(i,:)/sum(squareBipartiteGraphToVoteOn(i,:));
% $$$ end
