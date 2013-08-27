numCorrect = [15 25 31 37 41 48 53 53];
numIncorrect = [0 0 1 1 2 6 9 13];
plot(1:length(numCorrect),numCorrect, 'b+-', 1:length(numCorrect), numIncorrect, 'ro-',1:length(numCorrect),100*numCorrect./(numCorrect+numIncorrect), 'gx-',1:length(numCorrect),(numCorrect+numIncorrect)*100/70, 'k^-')
legend('numCorrect','numIncorrect','Assignment Accuracy','Percent assigned')
xlabel('iteration index')
%ylabel('number of assignments');
hold on
