%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Yasmin Fathy (Fathy.Yasmin@gmail.com)
% Last Updated: 4 June, 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Alg5_detection_increasing
% The new algorithm: this version the window size of the slow filter is
% a fixed window
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc       
clear      
close all  

%% Input Parameters to generate Input (signal) channel
channels_num = 1;
segments_num = 10;
variance = 1;
correlation_value = 0;
plot_flag = 0; % don't show/save figures of generated signal

learning_rate = 0.1; 

% window size for slow window (slow filter)
w_s = 50;  

% window size for fast window (fast filter)
w_f = 4; 

% detection threshold to check if hypothesis exceeds threshold, 
% a change is detected
threshold_h= 0.8;

% threshold_detection: a threshold to consider a change as a true positive 
threshold_detection = 50;

iter =1;
iterations = 5000;
tic % to calculate the execution time in seconds
for iter = 1:iterations
    iter

% x : input single-channel signal 
[x, corr_x, time_detection]= randomPieceWiseGenerator_fast(channels_num, segments_num, variance, correlation_value, plot_flag);


% moving average for fast filter
moving_avg_f = zeros(length(x),1); 

% moving average for slow filter
moving_avg_s = zeros(length(x),1);

lambda = 0; 
y = zeros(length(x),1); % estiamted signal
e = zeros(length(x),1); % error

% Calculate the fast filter
for k = w_f: length(x)
    moving_avg_f(k) = mean(x(k-w_f+1:k)); % moving average fixed window
end

% start index for slow filter
start = 1;

% k is the index over the signal x for slow filter
k = start + w_s ; 

% detection values = 1 where a change is detected and o otherwise
detection = zeros(length(x), 1);
while k < length(x) - 1    
    % calculate moving average for slow filter with an increasing widow
    moving_avg_s(k) = mean(x(start:k));% moving average increasing window
    
    % estimated signal
    y(k) = (lambda * moving_avg_f(k)) + ((1-lambda) * moving_avg_s(k));
    % weight is the diff between the two moving average
    w = moving_avg_f(k) -  moving_avg_s(k);
      
    e(k) = x(k+1) - y(k);
    % estimation of the mixing parameter of two filters in on-line fashion
    lambda = lambda + (learning_rate * e(k) * w);
    h(k) = lambda;
    if lambda > 1
       lambda = 1;      
    end
    
    if lambda < 0
        lambda = 0;
    end
    h(k) = lambda;
    % A decision rule: threshold for detetcing if change has occurred
    if lambda > threshold_h  && k<length(x)-1-w_s
        start = k;
        detection(k) = 1;
        % reset the starting-point (i.e. forgetting the factor)
        k = k + w_s;
        lambda = 0;
    else
         k = k + 1;       
    end
    
end

% get where changes are detected
nd_idx = find(detection==1);

% calculate mean square error
MSE_e(iter) = mean(e.^2); 

% counter
false_positive = 0;

% to calculate the delay if exists between a change 
delay= NaN; 
idx=1; % index for delay matrix "delay"

% to save the true positive-change detection time (i.e. the 
detection_tp = zeros(1,length(time_detection));


% it includes the elements in time_detection that have been detected
Detection=zeros(1,length(time_detection));

for i= 1:length(nd_idx)
    % flag to indicate if the current nc_idx(i) is a true positive
    is_tp = 0; 
    for j = 1: length(time_detection)
        % check if a change detected by algo (is in nd_idx) matches one of
        % changes in time_detection 
        if nd_idx(i) >= time_detection(j) && nd_idx(i) <=(time_detection(j) + threshold_detection)
           
           % match the current nd_idx(i) to the nearest value of 
           % time_detection
           [~,I] = min(abs(nd_idx(i)-time_detection));
           
           % mark which time_detection index matches the current nd_idx
           Detection(I) = 1;
                      
           % detect how much delay between the actucal change 
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
if isnan(delay)
    delay_avg = 0;
else
    delay_avg(iter) = mean(delay);
end
true_positive = sum(Detection);
x_len = length(x); % length of (generated) Input signal


% ERROR II type: FNR
false_negative = length(time_detection) - true_positive;
false_negative_rate(iter) = (false_negative)/(true_positive + false_negative);

true_positive_rate(iter)  = true_positive/(true_positive + false_negative);


% true negative rate is Specifivity
true_negative = x_len - true_positive;
true_negative_rate(iter) = true_negative/(true_negative + false_positive);


%ERROR I type: FPR
false_positive_rate(iter) = false_positive/(false_positive + true_negative);

end

FPR = mean(false_positive_rate) % should be minimal
FNR = mean(false_negative_rate) % should be minimal
delayR = mean(delay_avg) % should be minimal

timeElapsed = toc % to calculate execution time in seconds

%% plot all results
fig = figure;
plot(x),xlim([0,length(x)]),grid on;
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
title('Proposed approach: single sensor with an increasing window');
xlabel('n'),ylabel('x[n]');
legend('Location','southoutside')
legend([h1, h2], {'Actual', 'TP'});

% add title to the figure window
set(gcf,'name','Proposed approach with an increasing window','numbertitle','off');
savefig(fig,'single_proposed_increasing_2d.fig');

save('single_proposed_increasing_2d.mat','iterations','channels_num','segments_num','variance','correlation_value',...
'h','threshold_detection','x','time_detection','Detection','false_negative',...
'false_positive','true_positive','true_negative','delay','false_positive_rate',...
    'false_negative_rate','delayR','FDR','false_discovery_rate','FPR','FNR');
