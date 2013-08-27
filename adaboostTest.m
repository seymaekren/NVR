function adaboostTest

dbstop if error;
dbstop if warning;
    
%[TrainData, TrainLabels]     =
%convertFromLibSvmFormat('balancedData_scaled.txt');
%[TrainData, TrainLabels]     =
%convertFromLibSvmFormat('half_data.txt.scale');
[TrainData, TrainLabels]       = convertFromLibSvmFormat('data_scaled.txt');

[ControlData1, ControlLabels1] = convertFromLibSvmFormat('filteredCorrectScores.txt.scale');
[ControlData2, ControlLabels2] = ...
    convertFromLibSvmFormat('filteredIncorrectScores.txt.scale');

MaxIter = 100; % boosting iterations

% Step3: constructing weak learner
weak_learner = tree_node_w(3); % pass the number of tree splits to the constructor

% Step4: training with Gentle AdaBoost
[RLearners RWeights] = GentleAdaBoost(weak_learner, TrainData, TrainLabels, MaxIter);

% Step5: training with Modest AdaBoost
[MLearners MWeights] = ModestAdaBoost(weak_learner, TrainData, TrainLabels, MaxIter);

ControlData          = ControlData1; ControlLabels = ControlLabels1;

% Step6: evaluating on control set
ResultR = sign(Classify(RLearners, RWeights, ControlData));

Score_CorrectAssignments = Classify(RLearners, RWeights, ControlData);

ResultM = sign(Classify(MLearners, MWeights, ControlData));

% Step7: calculating error
ErrorR  = sum(ControlLabels ~= ResultR) / length(ControlLabels)

ErrorM  = sum(ControlLabels ~= ResultM) / length(ControlLabels)

ControlData = ControlData2; ControlLabels = ControlLabels2;

% Step6: evaluating on control set
ResultR = sign(Classify(RLearners, RWeights, ControlData));

Score_IncorrectAssignments = Classify(RLearners, RWeights, ControlData);

ResultM = sign(Classify(MLearners, MWeights, ControlData));

% Step7: calculating error
ErrorR  = sum(ControlLabels ~= ResultR) / length(ControlLabels)

ErrorM  = sum(ControlLabels ~= ResultM) / length(ControlLabels)


fid = fopen('adaboostScore_full_CorrectAssignments.txt', 'w');
for i = 1:length(Score_CorrectAssignments)
  fprintf(fid, '%f\n', Score_CorrectAssignments(i));
end
fclose(fid);

fid = fopen('adaboostScore_full_IncorrectAssignments.txt', 'w');
for i = 1:length(Score_IncorrectAssignments)
  fprintf(fid, '%f\n', Score_IncorrectAssignments(i));
end
fclose(fid);

keyboard

function [data, labels]     = convertFromLibSvmFormat(dataFilename);

%this is to compute the weights from the SVM model file.
[alpha xi1s xi2s xi3s xi4s] = textread(dataFilename,'%f %s %s %s %s');

data   = zeros(4,length(alpha));
labels = zeros(1, length(alpha));

for i = 1:length(alpha)

  [T,R] = strtok(xi1s(i),':');
  xi1   = strtok(R,':');
  xi1   = str2num(xi1{1});
  
  [T,R] = strtok(xi2s(i),':');
  xi2   = strtok(R,':');
  xi2   = str2num(xi2{1});
  
  
  [T,R] = strtok(xi3s(i),':');
  xi3   = strtok(R,':');
  xi3   = str2num(xi3{1});
  
  
  [T,R] = strtok(xi4s(i),':');
  xi4   = strtok(R,':');
  xi4   = str2num(xi4{1});
  
  
  xi        = [xi1 xi2 xi3 xi4];
  data(:,i) = xi';
  labels(i) = alpha(i);
end