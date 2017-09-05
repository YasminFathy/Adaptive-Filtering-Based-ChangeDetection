
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Yasmin Fathy (Fathy.Yasmin@gmail.com)
% Last Updated: 5 Sept, 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Proposed Algorithm: An On-line Adaptive Algorithm for Change
%                     Detection in Streaming Sensory Data
% The new algorithm: this version the window size of the slow filter is
% an increasing window in distributed, cooperative manner between the signals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
clc        % clear command window
clear      % clear all variables
close all  % close all figures

%% Input Parameters to generate Input (signal) channel
sensors_num = 10;

channels_num = sensors_num;

segments_num = 10;
variance = 1;
correlation_value = 0;
plot_flag = 0; % don't how/save figures of generated signal

%% Input parameters for Alg5_detection_fixed
% learning rate for the weights
learning_rate = 0.1; 

% window size for slow window (slow filter)
w_s = 50;

% window size for fast window (fast filter)
w_f = 4; %10; %5; 

% detection threshold to check if hypothesis exceeds threshold, 
% a change is detected
threshold_h= 0.6; % this value affects the delay of detecting a change

% threshold_detection: a threshold to consider a change as a true positive 
% if the change is detected within a margin of that threshold; this for
% evaluation not for the alg itself
% threshold_detection = w_s;
threshold_detection = 50;

iter =1;
iterations = 1;
tic % to calculate the execution time in seconds
for iter = 1:iterations
    iter

% x : input single-channel signal 
% time_detection : the actual changes in the input generated signal
% corr_x : correlated multi-channel signal
% time_detection : time where changes occured
% corr_mat: correlation matrix between different channels
% mu : list of mu of the segments
% length_segments: list of the length of each segment

[x, corr_x, time_detection]= randomPieceWiseGenerator(channels_num, segments_num, variance, correlation_value, plot_flag);
if correlation_value > 0
    x = corr_x; 
end       
    % Calculate the fast filter for all channels(sensors)
    % We assume that the fast filter values are calculated on parallel with
    % the slow filter
    for k = w_f: length(x)-1
        % just note:  the filter is online, so the value of moving_avg_f(k) is
        % calculated based on the values we have seen it so far until the
        % current t
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
    for sn = 1:sensors_num
        % calculate moving average for slow filter with an increasing widow
        % just note:  the filter is online, so the value of moving_avg_s(k) is
        % calculated based on the values we have seen it so far until the
        % current t
        moving_avg_s(sn,k) = mean(x(sn,start:k));% moving average increasing window
        % estimated signal
        % output of overall filters (fast and slow filters)
        % y(k) is the total weight
        y(sn) = (lambda(sn) * moving_avg_f(sn,k)) + ((1-lambda(sn)) * moving_avg_s(sn,k));  
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
    % A decision rule: threshold for detetcing if change has occurred
    if common_lambda > threshold_h  && k<length(x)-(1-w_s)
        start = k;
        detection(k) = 1;
        % reset the starting-point (i.e. forgetting the factor)
        k = k + w_s;
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

% counter
false_positive = 0;

% to calculate the delay if exists between a change 
% (that is actually happened (i.e. is in the generated signal) which is in 
% time_detection) and the same one if it is detected by the algorithm 
% (i.e. is in nd_idx)
delay= NaN; 
idx=1; % index for delay matrix "delay"

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

end

FPR = mean(false_positive_rate) % should be minimal
FNR = mean(false_negative_rate) % should be minimal
delayR = mean(delay_avg) % should be minimal
MSE = mean(MSE)  % should be minimal
TPR = mean(true_positive_rate)
TNR = mean(true_negative_rate)

timeElapsed = toc % to calculate execution time in seconds

%% plot all results
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
%title('Abrupt Detection in Signal x[n]');
title('Proposed approach: multiple sensors with an increasing window');
xlabel('n'),ylabel('x[n]');
legend('Location','southoutside')
legend([h1, h2], {'Actual', 'TP'});

% add title to the figure window
set(gcf,'name','Proposed approach with an increasing window','numbertitle','off');

fname1=sprintf('multi_proposed_increasing_2d_%dchs_%dcorr_500LR_%d_h%d.fig',sensors_num,correlation_value,learning_rate,threshold_h);
%'multi_proposed_increasing_10chs_0_1corr.fig'
savefig(fig,fname1);

fname2=sprintf('multi_proposed_increasing_2d_%dchs_%dcorr_500LR_%d_h%d.mat',sensors_num,correlation_value,learning_rate,threshold_h);
% save experiments results and set-up in mat file
save(fname2,'iterations','channels_num','segments_num','variance','correlation_value',...
'h','threshold_detection','x','time_detection','Detection','false_negative',...
'false_positive','true_positive','true_negative','delay','true_positive_rate','true_negative_rate','false_positive_rate',...
    'false_negative_rate','delayR',,'FPR','FNR',...
    'TNR','TPR','w_s','w_f','threshold_h','learning_rate','threshold_detection');
