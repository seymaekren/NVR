function printAssignments (filename, assignments)

fout          = fopen  (filename ,'w');
fprintf(fout, '%d %d\n', assignments');
fclose (fout);
