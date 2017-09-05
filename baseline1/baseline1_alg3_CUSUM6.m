%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Yasmin Fathy (Fathy.Yasmin@gmail.com)
% Last Updated: 25 July, 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Alg3_CUSUM6

% Algorithm 3: CUSUM with recursive form p. 6 (it is implemented it in
% a sub-optimal form from p.11 in same reference)
% ref: Granjon, P. (2012). The cusum algorithm a small review. Gipsa-Lab, 
% Grenoble, France, Team SAIGA.
% The implementaion is based on calaulating instantaneous log-likelihood
% s(n) using equ 1.19 in p.10 in the same reference.


clc        
clear      
close all  

channels_num = 1;
segments_num = 10;
variance = 1;
correlation_value = 0;
plot_flag = 0; % don't show/save figures of generated signal
h = 8; 
d = 2; 

threshold_detection = 50;

iter =1;
iterations = 1;
tic % to calculate the execution time in seconds
for iter = 1:iterations
    iter
[x, corr_x, time_detection]= randomPieceWiseGenerator(channels_num, segments_num, variance, correlation_value, plot_flag);

nc = zeros(length(x),1); 
nd = zeros(length(x),1); 


s = zeros(length(x), 1); % instantaneous log likelihood
S = zeros(length(x), 1); % cumulative sum
G = zeros(length(x), 1); % decision function

start = 2;

for k=start:length(x)
    mu = mean(x(start-1:k));  % current mean
    sigma= std(x(start-1:k)); % current std
    % calculate instantaneous log-likelihood
    s(k) = (d/sigma^2) * (x(k) - mu- (d/2));
    % calcualte cumulative sum
    S(k)= S(k-1) + s(k);
    % calculate the decision function
    G(k) = max(G(k-1) + s(k), 0); % supremum
    
    % compare the decision function with the threshold
    if G(k) > h % H1 occurs, H0 not
        nd(k) = 1; % detection happens after the changes have happened
        
        % estimate change location
        [~,I] = min(S(1:k-1)); % I is the indices where it has the min
        nc(I) = 1;        
        
        % Reset the algorithm
        start = k;
        S(start-1) = 0;
        G(start-1) = 0;        
    end
end
% get where changes are detected
nc_idx=find(nc==1); % estimation change location

% nd_idx: indexes of all changes detected by the algorithm 
% find is to get indices that have only 1 which t values.
nd_idx=find(nd==1); % detection of the changes



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

FPR = mean(false_positive_rate) % should be minimal
FNR = mean(false_negative_rate) % should be minimal
delayR = mean(delay) % should be minimal


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
title('One-sided CUSUM: single sensor');
xlabel('n'),ylabel('x[n]');
legend('Location','Southoutside')
legend([h1, h2], {'Actual', 'TP'});

% add title to the figure window
set(gcf,'name','Baseline 1','numbertitle','off');

savefig(fig,'baseline1_1d.fig');

% save experiments results and set-up in mat file
save('baseline1_1d.mat','iterations','channels_num','segments_num','variance','correlation_value',...
'h','d','threshold_detection','x','time_detection','Detection','false_negative',...
'false_positive','true_positive','true_negative','delay','false_positive_rate',...
    'false_negative_rate','delayR','FPR','FNR');
