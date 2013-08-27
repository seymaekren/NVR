function retval = MC_Accept(prevScore, score)

kT = 1;

score     = score;
prevScore = prevScore;

if (score < prevScore)
  retval = 1;
  return
end

prob = exp(-(score - prevScore)/kT);
if (rand() < prob)
  retval = 1;
else
  retval = 0;
end
 
%keyboard