load parsedBackboneNOEs.txt
parsedBackboneNOEs = parsedBackboneNOEs+720;
fout  = fopen('parsedBackboneNOEs_renumbered.txt','w');
for i = 1:size(parsedBackboneNOEs,1)
  fprintf(fout, '%d \t%d\n',parsedBackboneNOEs(i,1),parsedBackboneNOEs(i,2));
end
fclose(fout);

