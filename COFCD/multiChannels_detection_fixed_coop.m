%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Yasmin Fathy (Fathy.Yasmin@gmail.com)
% Last Updated: 5 Sept, 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Proposed Algorithm: An On-line Adaptive Algorithm for Change
%                     Detection in Streaming Sensory Data
% The new algorithm: this version the window size of the slow filter is
% a fixed window in distributed, cooperative manner between the signals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc        
clear     
close all  

%% Input Parameters to generate Input (signal) channel
sensors_num = 10;

channels_num = sensors_num;

segments_num = 10;
variance = 1;
correlation_value = 0.1; 
plot_flag = 0; % don't show/save figures of generated signal
learning_rate = 0.1; 

% window size for slow window (slow filter)
w_s = 50;  

% window size for fast window (fast filter)
w_f = 4; 
threshold_h= 0.6;
threshold_detection = 50;

iter =1;
iterations = 1000;
tic % to calculate the execution time in seconds
for iter = 1:iterations
    iter

[x, corr_x, time_detection]= randomPieceWiseGenerator(channels_num, segments_num, variance, correlation_value, plot_flag);
if correlation_value > 0
   x = corr_x; 

end  
    for k = w_f: length(x)-1
        moving_avg_f(:,k) = mean(x(:,k-w_f+1:k),2); % moving average fixed window
    end
    
   detection = zeros(length(x), 1);   
   lambda = zeros(sensors_num, 1); % weight for each sensor
   
    common_lambda = 0; % common lambda across channels/sensors
    start = 1; % start index for slow filter
    % k is the index over the signal x for slow filter
    k = start + w_s ; 
    
 while k < length(x) - 1      
    % iterate for each channel/sensor   
    for sn = 1:sensors_num         
        moving_avg_s(sn,k) = mean(x(sn,k-w_s+1:k)); % moving average with fixed window
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
        lambda = zeros(sensors_num, 1); % weight for each sensor
        common_lambda = 0;
    else
         k = k + 1;       
    end    
    % avg mean sqaure error
    mse(k) = mean(e);

 end
 nd_idx=find(detection==1); % detection of the changes
 
% calculate mean square error
MSE(iter) = mean(mse.^2); 


false_positive = 0;
delay= NaN; 
idx=1; % index for delay matrix "delay"

detection_tp = zeros(1,length(time_detection));
Detection=zeros(1,length(time_detection));

for i= 1:length(nd_idx)
    is_tp = 0; 
    for j = 1: length(time_detection)
        if nd_idx(i) >= time_detection(j) && nd_idx(i) <=(time_detection(j) + threshold_detection)
           [~,I] = min(abs(nd_idx(i)-time_detection));
           Detection(I) = 1;
           delay(idx) = abs(nd_idx(i) - time_detection(I));
           idx = idx + 1;
           detection_tp(I) = nd_idx(i);
           is_tp = 1;   
        end 
    end
    if is_tp == 0
          false_positive = false_positive + 1;
    end
end

if isnan(delay)
    delay_avg = 0;
else
    delay_avg(iter) = mean(delay);
end

true_positive = sum(Detection);
x_len = length(x); 

false_negative = length(time_detection) - true_positive;
false_negative_rate(iter) = (false_negative)/(true_positive + false_negative);
true_positive_rate(iter)  = true_positive/(true_positive + false_negative);
true_negative = x_len - true_positive;
true_negative_rate(iter) = true_negative/(true_negative + false_positive);

false_positive_rate(iter) = false_positive/(false_positive + true_negative);


end

FPR = mean(false_positive_rate) 
FNR = mean(false_negative_rate) 
delayR = mean(delay_avg) 
MSE = mean(MSE)  
TPR = mean(true_positive_rate)
TNR = mean(true_negative_rate)

timeElapsed = toc % to calculate execution time in seconds

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
fname1=sprintf('multi_proposed_fixed_1d_%dchs_%dcorr.fig',sensors_num,correlation_value);
%'multi_proposed_increasing_10chs_0_1corr.fig'

savefig(fig,fname1);
fname2=sprintf('multi_proposed_fixed_2d_%dchs_%dcorr.mat',sensors_num,correlation_value);

% save experiments results and set-up in mat file
save(fname2,'iterations','channels_num','segments_num','variance','correlation_value',...
'h','threshold_detection','x','time_detection','Detection','false_negative',...
'false_positive','true_positive','true_negative','delay','true_positive_rate','true_negative_rate','false_positive_rate',...
    'false_negative_rate','delayR','FPR','FNR',...
    'TNR','TPR','w_s','w_f','threshold_h','learning_rate','threshold_detection');
