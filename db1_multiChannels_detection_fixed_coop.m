%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Yasmin Fathy (Fathy.Yasmin@gmail.com)
% Last Updated: 25 July, 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Proposed Algorithm: An On-line Adaptive Algorithm for Change
%                     Detection in Streaming Sensory Data
% The new algorithm: this version the window size of the slow filter is
% a fixed window in distributed, cooperative manner between the signals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc        % clear command window
clear      % clear all variables
close all  % close all figures

dataset1;


%% Input parameters for Alg5_detection_fixed
% learning rate for the weights
learning_rate = 0.1; 

% window size for slow window (slow filter)
w_s = 40;  %40; %50

% window size for fast window (fast filter)
w_f = 5; %10; %5; 

% detection threshold to check if hypothesis exceeds threshold, 
% a change is detected
threshold_h= 0.6;

% threshold_detection: a threshold to consider a change as a true positive 
% if the change is detected within a margin of that threshold; this for
% evaluation not for the alg itself
% threshold_detection = w_s;
%threshold_detection = 50;
threshold_detection = 50;

iter =1;
iterations = 1;
tic % to calculate the execution time in seconds
for iter = 1:iterations
    iter


%     moving average for fast filter
%     moving_avg_f = zeros(sensors_num,length(x));     
    
%     moving average for slow filter
%     moving_avg_s = zeros(length(x),1);
    
%     y = zeros(length(x),1); % estiamted signal
%     e = zeros(length(x),1); % error
    
    % Calculate the fast filter for all channels(sensors)
    % We assume that the fast filter values are calculated on parallel with
    % the slow filter
    for k = w_f: length(x)-1
        % just note:  the filter is online, so the value of moving_avg_f(k) is
        % calculated based on the values we have seen it so far until the
        % current t
        % mean(A,2) is a column vector containing the mean of each row.
        moving_avg_f(:,k) = mean(x(:,k-w_f+1:k),2); % moving average fixed window
    end

    
   % detection values = 1 where a change is detected and o otherwise
   detection = zeros(length(x), 1);   
   lambda = zeros(sensors_num, 1); % weight for each sensor
   
    common_lambda = 0; % common lambda across channels/sensors
    start = 1; % start index for slow filter
    % k is the index over the signal x for slow filter
    k = start + w_s ; 
    
 while k < length(x) - 1      
    % iterate for each channel/sensor 
    %for sensor = 1:sensors_num
    %common_lambda = mean(lambda);
    %h(k) = common_lambda;
    
    for sn = 1:sensors_num         
        % calculate moving average for slow filter with a fixed widow
        % just note:  the filter is online, so the value of moving_avg_s(k) is
        % calculated based on the values we have seen it so far until the
        % current t
        moving_avg_s(sn,k) = mean(x(sn,k-w_s+1:k)); % moving average with fixed window

         % estimated signal
        % output of overall filters (fast and slow filters)
        % y(k) is the total weight
        y(sn) = (lambda(sn) * moving_avg_f(sn,k)) + ((1-lambda(sn)) * moving_avg_s(sn,k));
        %e(k) = x(k+1) - y(k);    
        e(sn) = x(sn,k+1) - y(sn);
        % weight is the diff between the two moving average
        w = moving_avg_f(sn,k) -  moving_avg_s(sn,k);
        % estimation of the mixing parameter of two filters in on-line fashion
        lambda(sn) = common_lambda + (learning_rate * e(sn) * w);
        % check on common(accumulated) lambda
        if lambda(sn) > 1
               lambda(sn) = 1;      
        end

        if lambda(sn) < 0
            lambda(sn) = 0;
        end
    end
    common_lambda = mean(lambda);
    h(k) = common_lambda;

    %h(k) = common_lambda;
    % A decision rule: threshold for detetcing if change has occurred
    if common_lambda > threshold_h  && k<length(x)-(1-w_s)
        start = k;
        detection(k) = 1;
        % reset the starting-point (i.e. forgetting the factor)
        k = k + w_s;
        %h(k) = 0;
        %lambda = 0;
        lambda = zeros(sensors_num, 1); % weight for each sensor
        common_lambda = 0;
    else
         k = k + 1;       
    end    
    % avg mean sqaure error
    mse(k) = mean(e);

 end
 % get where changes are detected
 % nd_idx: indexes of all changes detected by the algorithm 
 % find is to get indices that have only 1 which t values.
 nd_idx=find(detection==1); % detection of the changes
 
% calculate mean square error
MSE(iter) = mean(mse.^2); 
%%get overlapped indxs
 
% % get lengths of all nd_idx for all channels 
% length_ch= cellfun(@length,nd_idx_arr);
% % to have chanels with same length, we get the min length of all channels
% min_length = min(length_ch);
% % % intialise all indices by zeros
% nd_idx_overlapped= zeros(min_length,length(length_ch));   
% % add 10 to the values on each channel (each cell in nd_idx_distributed1)
% for j=1: length(length_ch)
%        nd_idx_overlapped(:,j) = nd_idx_arr{j}(1:min_length); 
% 
% end

% %get lengths of all nd_idx for all channels in nd_idx3
% length_ch= cellfun(@length,nd_idx_arr);
% % to have chanels with same length, we get the min length of all channels
% %min_length = min(length_ch);
% %min_length = mode(length_ch);
% min_length = median(length_ch);
% % intialise all indices by zeros
% nd_idx_overlapped= zeros(min_length,length(length_ch));    
% % add 50 to the values on each channel (each cell in nd_idx_distributed1)
% for j=1: length(length_ch)
%     % in case the current length < the agreed length
%     ch_length = length(nd_idx_arr{j});
%     if ch_length < min_length
%         for k = (ch_length+1) : min_length
%             nd_idx_arr{j}(k) = 0;
%         end 
%     end
%    nd_idx_overlapped(:,j) = nd_idx_arr{j}(1:min_length); 
% end
% 
% 
% % get the overlapped value of the changes across each row in the matrix
% %nd_idx= mode(nd_idx_overlapped,2);
% nd_idx= median(nd_idx_overlapped,2);
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
% mean(delay_avg) % should be minimal
% mean(MSE)  % should be minimal

FPR = mean(false_positive_rate) % should be minimal
FNR = mean(false_negative_rate) % should be minimal
delayR = mean(delay_avg) % should be minimal
MSE = mean(MSE)  % should be minimal
FDR = nanmean(false_discovery_rate)
TPR = mean(true_positive_rate)
TNR = mean(true_negative_rate)

timeElapsed = toc % to calculate execution time in seconds

%mean(true_positive_rate); % should be maximised
%mean(true_negative_rate); % should be maximised
%mean(precision);          % should be maximised
%mean(false_discovery_rate); % should be minimal

%% plot all results
%plot(detection);
fig = figure;
plot(x'),xlim([0,length(x)]),grid on;
ylim = get(gca,'YLim');
hold on

% plot where changes are actually happened from the generated signal
for i=1:length(time_detection)
    h1=line([time_detection(i) time_detection(i)], ylim,'Color', [0 0 0]);
end

% detection_tp includes all
detection_tp = detection_tp(detection_tp>0);

% plot where changes are detected by the algorithm
for i=1:length(detection_tp)
    h2=line([detection_tp(i) detection_tp(i)], ylim, 'LineStyle','--','Color', [1 0 0]);
end
hold off

title('Proposed approach: multiple sensors with a fixed window');
xlabel('n'),ylabel('x[n]');
legend('Location','southoutside')
legend([h1, h2], {'Actual', 'TP'});

% add title to the figure window
set(gcf,'name','Proposed approach with a fixed window','numbertitle','off');

% fname1=sprintf('realwordldata_proposed_fixed_%dchs_LR_%d_h%d.fig',sensors_num,learning_rate,threshold_h);
fname1=sprintf('realwordldata_db1_proposed_fixed.fig');
%'multi_proposed_increasing_10chs_0_1corr.fig'

savefig(fig,fname1);
