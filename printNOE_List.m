function printNOE_List(NOES)

fid = fopen('NOE_List.txt','w');
fprintf(1, 'printing to NOE_List.txt\n');

filename2 = 'rawNOE_List_peakIndices.txt';
fprintf(1, 'printing to %s\n',filename2);
fid2 = fopen(filename2,'w');


for peak1Index = 1:size(NOES,1)
  listOfNOEs = find(NOES(peak1Index,:));
  if (~isempty(listOfNOEs))
    printThisPeak             = 0;
    listOfPeaksToBePrinted    = [];
    for peak2Index = 1:size(NOES,2)
      %	fprintf(fid, '%d ', NOES(ROWIN(rowIndex),
      %	COLIN(columnIndex)));
      if (peak2Index <= peak1Index)
	continue;
      end
      if (NOES(peak1Index,peak2Index) == 1)
	if (printThisPeak == 0)
%	  fprintf(fid, '%d %d ', length(listOfNOEs), peak1Index);
	  listOfPeaksToBePrinted = [peak1Index];
	  printThisPeak          = 1;
	end
	listOfPeaksToBePrinted = [listOfPeaksToBePrinted peak2Index];
%	fprintf(fid, ' %d ',peak2Index);
        fprintf(fid2, '%d %d\n', peak1Index, peak2Index);
      end
    end

    if (printThisPeak)
      assert (length(listOfPeaksToBePrinted)>1);
      fprintf(fid, '%d %d ', length(listOfPeaksToBePrinted)-1, listOfPeaksToBePrinted(1));
      for peaksToBePrintedIndex = 2:length(listOfPeaksToBePrinted)
	fprintf(fid, '%d ', listOfPeaksToBePrinted(peaksToBePrintedIndex));
      end
      fprintf(fid, '\n');
    end
  end
end

fclose(fid);
fclose(fid2);
