load chRdc.txt
chRdc(:,1) = chRdc(:,1)+720;
fout  = fopen('chRdc_renumbered.txt','w');
for i = 1:size(chRdc,1)
  fprintf(fout, '%d \t%f \t%f \t%f\n',chRdc(i,1),chRdc(i,2),chRdc(i,3),chRdc(i,4));
end

load nhRdc.txt
nhRdc(:,1) = nhRdc(:,1)+720;
fout  = fopen('nhRdc_renumbered.txt','w');
for i = 1:size(nhRdc,1)
  fprintf(fout, '%d \t%f \t%f \t%f\n',nhRdc(i,1),nhRdc(i,2),nhRdc(i,3),nhRdc(i,4));
end
fclose(fout);

load parsedResonances.txt
parsedResonances(:,1)  = parsedResonances(:,1)+720;
fout  = fopen('parsedResonances_renumbered.txt','w');
for i = 1:size(parsedResonances,1)
  fprintf(fout, '%d %.2f %.2f\n', parsedResonances(i,1),parsedResonances(i,2),parsedResonances(i,3));
end
fclose(fout);