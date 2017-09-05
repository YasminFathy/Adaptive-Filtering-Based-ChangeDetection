%clear;
%%

%%%% Added by Yasmin to have a fair comparion;using the same 
% generated function


iter =1;
iterations = 1;
threshold_detection = 50;

%  f = fullfile('../../real-world-dataset', 'database1.m');
 run('real-world-dataset/dataset1.m');

tic % to calculate the execution time in seconds
for iter = 1:iterations
    iter

    y = x;

%%%% Finish updated by Yasmin

%load logwell.mat

% choice of alpha
alpha = 0.01;%.01;

n = 50;  %50% no of samles
k = 10; %50;  %10; % window size
% subplot(2,1,1)
% plot(y')
% axis([1,size(y,2),-inf,inf])

%%
% Forward detectino
score1 = change_detection(y,n,k,alpha,5);
% Backword detection
score2 = change_detection(y(:,end:-1:1),n,k,alpha,5);


% Since the alg plots the pdf of the changes, to have a fair comparions
% with our algorithm, I'd like to get first few highest detections
data= [score1 + score2];
[peaks_values, peaks_locs] = findpeaks(data);
% Sort them in a descending order to find the largest peaks.
[sorted_values, sorted_inds] = sort(peaks_values, 'descend');
% get original locations of the peaks
originalLocations = peaks_locs(sorted_inds);
% get the highest n peaks where n = the number of actual changes generated 
% in the data
%detected_locations = originalLocations(1:(length(time_detection)));
% find cloest values from detected_locations to each value in
% time_detection

for i=1:length(time_detection)
    value = time_detection(i);
    tmp = abs(originalLocations-value);
    [val idx] = min(tmp); %index of closest value
    closest = originalLocations(idx); %closest value
    det_locations(i) = closest;
    
end
%nd_idx = det_locations;
nd_idx = originalLocations;
%% evaluation
% calculate the false positive, false negative and delay
%%% 
% A false positive is a result that indicates a given condition has been 
% fulfilled, when it has not. 
% A change is detected, however it did not happen (it was not exist in
% the generated signal)

% Check:
%http://web.stanford.edu/~rjohari/teaching/notes/226_lecture8_prediction.pdf

% A false negative is where a test result indicates that a condition 
% failed, while it was successful
% A change was not detected, however it did happen (it was exist in the 
% generated signal)
%
% A true positive: is where a change happened (in the generated signal) 
% and it was detected by the algorithm
%
% A true negative: is where a change did not happen (in the generated 
% signal) and it was not detected by the algorithm
%%%


% counter
false_positive = 0;

% to calculate the delay if exists between a change 
% (that is actually happened (i.e. is in the generated signal) which is in 
% time_detection) and the same one if it is detected by the algorithm 
% (i.e. is in nd_idx)
delay= NaN; 
idx=1; % index for delay matrix "delay"

% to make each nd_idx to only map to one element in time_detection
% Note: this is because often once the algorithm detects a change, it keeps
% detecting it for some next following signals/samples, so I have to
% detect that change only once and it should just match one element in
% time_detection

%dt=diff(nd_idx)>5;
%nd_idx=nd_idx(dt,:);

% to save the true positive-change detection time (i.e. the 
% changes detected correctly by the algorithm) and match 
% elements in time_detection
detection_tp = zeros(1,length(time_detection));


% it includes the elements in time_detection that have been detected
% by nd_idx (output of the alg), it will have 1 if it's detected 
% and 0 if it is not
Detection=zeros(1,length(time_detection));

for i= 1:length(nd_idx)
    % flag to indicate if the current nc_idx(i) is a true positive
    % it should be reinitialised within the loop to indicate if the current 
    % nc_idx1(i) is true positive or not.
    is_tp = 0; 
    for j = 1: length(time_detection)
        % check if a change detected by algo (is in nd_idx) matches one of
        % changes in time_detection (which are from generated signal)
        % and detect the delay if exists
        if nd_idx(i) >= time_detection(j) && nd_idx(i) <=(time_detection(j) + threshold_detection)
           
           % match the current nd_idx(i) to the nearest value of 
           % time_detection (i.e. with min distance)
           % Note: I here is index of one of elements in time_detection
           [~,I] = min(abs(nd_idx(i)-time_detection));
           
           % mark which time_detection index matches the current nd_idx
           % (true positive)
           Detection(I) = 1;
                      
           % detect how much delay between the actucal change 
           % time_detection(I) and the detected change by the algorithm
           % nc_idx(i)
           %delay(idx) = abs(nd_idx(i) - time_detection(j));
           delay(idx) = abs(nd_idx(i) - time_detection(I));
           
           % increase delay index (idx)
           idx = idx + 1;
           
           % to save which nd_idx matches the time_detection
           %detection_tp(idx) = nd_idx(i);
           detection_tp(I) = nd_idx(i);
  
           % flag that the current value (nc_idx(i)) is a true positive
           is_tp = 1;   
        end 
    end
    % a change is detected by the algorithm, but it did not really map to
    % any elements in time_detection
    if is_tp == 0
          false_positive = false_positive + 1;
    end
end

% get the average delay between the detected changes by the algorithm
% (i.e. are in nd_idx) and between the actual changes in the generated 
% signal (are in time_detection)
if isnan(delay)
    delay_avg = 0;
else
    delay_avg(iter) = mean(delay);
end
% calculate how many (actual changes)-elements in time_detection have  been 
% detected by the algorithm
true_positive = sum(Detection);
x_len = length(x); % length of (generated) Input signal


% ERROR II type: FNR
% false negative = #actaul_changes - #detected changes
%false_negative = abs(length(time_detection) - true_positive);
false_negative = length(time_detection) - true_positive;
% false_negative_rate = FN / (FN + TP)
false_negative_rate(iter) = (false_negative)/(true_positive + false_negative);

% true positive rate is Detection Rate/Sensitivity
% true_positive_rate = TP / (TP + FN)
%true_positive_rate  = sum(Detection)/(sum(Detection) + false_negative);
true_positive_rate(iter)  = true_positive/(true_positive + false_negative);


% true negative rate is Specifivity
% true_negative_rate = TN / (TN + FP)
% true negative is other values in x that are not true_positive
% true_positive which is the count of true positive that been detected 
% as mentioned above
true_negative = x_len - true_positive;
true_negative_rate(iter) = true_negative/(true_negative + false_positive);


%ERROR I type: FPR
% false positive rate is False Alarm Rate
% False Alarm Rate = FP / (FP + TN)
false_positive_rate(iter) = false_positive/(false_positive + true_negative);
%false_positive_rate = false_positive/(x_len);


% Precision = TP/(TP + FP)
precision(iter) = true_positive / (true_positive + false_positive);

% False discovery rate = FP / (TP+FP)
false_discovery_rate(iter) = false_positive/(true_positive + false_positive);

end

% mean(false_positive_rate) % should be minimal
% mean(false_negative_rate) % should be minimal
% mean(delay) % should be minimal

FPR = mean(false_positive_rate) % should be minimal
FNR = mean(false_negative_rate) % should be minimal
delayR = mean(delay) % should be minimal
FDR = nanmean(false_discovery_rate) % false discovery rate


timeElapsed = toc % to calculate execution time in seconds

%mean(true_positive_rate); % should be maximised
%mean(true_negative_rate); % should be maximised
%mean(precision);          % should be maximised
%mean(false_discovery_rate); % should be minimal


%% plots figure

fig = figure;
subplot(2,1,1);
%grid on;
%plot(y', 'linewidth',2);
plot(y');
%axis([-inf,size(y,2),-inf,inf])
title('Original Signal')
% Added by Yasmin: just add to plot actual changes in the data
% plot where changes are actually happened from the generated signal
hold on
ylim1 = get(gca,'YLim');
for i=1:length(time_detection)
    h1=line([time_detection(i) time_detection(i)], ylim1,'LineStyle','--','Color', [1 0 0]);
end
axis([-inf,size(y,2),-inf,inf])
grid on;
%%

subplot(2,1,2);
score2 = score2(end:-1:1);


% 2*n-2+k is the size of buffer zone
plot([zeros(1,2*n-2+k),score1 + score2], 'r-', 'linewidth',1,'Color', [0 0 0]);
title('Change-Point Scores')

% forward detection
%plot([zeros(1,2*n-2+k),score1], 'r-', 'linewidth',1,'Color', [0 0 1]);
% backword detection
%plot([zeros(1,2*n-2+k),score2], 'r-', 'linewidth',1,'Color', [1 0 0]);
%axis([-inf,size(y,2),-inf,inf])



% just plot points detected by algorithm
hold on;
ylim = get(gca,'YLim');

% for i=1:length(det_locations)
%     h1=line([det_locations(i) det_locations(i)], [0,ylim(2)],'linewidth',2, 'Color', [0 0 1]);
% end

for i=1:length(detection_tp)
    h1=line([detection_tp(i) detection_tp(i)], [0,ylim(2)],'linewidth',2, 'Color', [0 0 1]);
end

% plot where changes are actually happened from the generated signal
for i=1:length(time_detection)
    h1=line([time_detection(i) time_detection(i)], [0,ylim(2)], 'linewidth',2, 'LineStyle','--','Color', [1 0 0]);
end

axis([-inf,size(y,2),0,inf])
grid on;
savefig(fig,'realworld_db1_baseline3.fig');



