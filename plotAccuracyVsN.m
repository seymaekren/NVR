load Results_1UBI_sparse.txt
load Results_ff2.txt
load Results_hSRI_sparse.txt
figure;
numPeaks_UBI  = Results_1UBI_sparse(size(Results_1UBI_sparse,1),1);
numPeaks_ff2  = Results_ff2 (size(Results_ff2,1) ,1);
numPeaks_hSRI = Results_hSRI_sparse(size(Results_hSRI_sparse,1),1);

plot(100*Results_1UBI_sparse(:,1)/numPeaks_UBI, Results_1UBI_sparse(:,2), '*-', ...
     100*Results_ff2(:,1)/numPeaks_ff2, Results_ff2(:,2), 'o-', 100*Results_hSRI_sparse(:,1)/numPeaks_hSRI, Results_hSRI_sparse(:,2), 'x-');
legend('1UBI', 'ff2', 'hSRI');

title('Assignment Accuracy Vs. N');
xlabel('Percent of Assigned Peaks')
ylabel('Assignment Accuracy')
print -depsc AA_vs_N.eps
print -djpeg AA_vs_N.jpg