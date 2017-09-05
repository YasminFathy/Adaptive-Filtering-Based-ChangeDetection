%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Yasmin Fathy (Fathy.Yasmin@gmail.com)
% Last Updated: 25 July, 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Alg5_CUSUM17_increasing
% Algorithm 5: CUSUM with recursive form p. 17
% Version with an increasing window
% ref: Granjon, P. (2012). The cusum algorithm a small review. Gipsa-Lab, 
% Grenoble, France, Team SAIGA.
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

h = 8; 
d = 2; 
threshold_detection = 50;

iter =1;
iterations = 1000;
tic % to calculate the execution time in seconds
for iter = 1:iterations
    iter

[x, corr_x, time_detection]= randomPieceWiseGenerator(channels_num, segments_num, variance, correlation_value, plot_flag);
y=x;
% to detect increasing in a parameter
S_i = zeros(length(y),1);
s_i = zeros(length(y),1);
G_i = zeros(length(y),1);

% to detect decreasing in a parameter
S_d = zeros(length(y),1);
s_d = zeros(length(y),1);
G_d = zeros(length(y),1);

nc = zeros(length(y),1);
nd = zeros(length(y),1);

% intialisation
mu = zeros(length(y),1);
sigma = zeros(length(y),1);

S_i(1) = 0;
G_i(1) = 0;
S_d(1) = 0;
G_d(1) = 0;


start = 2;
for k=start:length(y)
    % calculate the estimation mu and sigma
    mu(k)= mean(y(start-1:k));
    sigma(k) = std(y(start-1:k))^2;
    
    % instantaneous log likelihood for increasing _i and decreasing _d
    s_i(k) = d/sigma(k) * (y(k) - mu(k) - (d/2));
    s_d(k) = -1 * d/sigma(k) * (y(k) - mu(k) + (d/2));
    
    % cumulative sum for increasing _i and decreasing _d
    S_i(k) = S_i(k-1) + s_i(k);
    S_d(k) = S_d(k-1) + s_d(k);
    
    % decision function for increasing _i and decreasing _d
    G_i(k) = max(G_i(k-1) + s_i(k),0);
    G_d(k) = max(G_d(k-1) + s_d(k),0);
    
    %G_i(k) = G_i(k-1) + s_i(k);
    %G_d(k) = G_d(k-1) + s_d(k);
    
    if G_i(k) > h || G_d(k) > h
        % detection happens after the changes have happened
        % n_d(k) = k;
        nd(k) = 1; 
        
        if G_i(k) > h
            % estimate change location in case of increasing
            nc= min(S_i(1:k-1));
            [~,I] = min(S_i(1:k-1)); % I is the indices where it has the min
            nc(I) = 1;        
            
        end
        if G_d(k) > h
            % estimate change location in case of decreasing
            nc = min(S_d(1:k-1));
            [~,D] = min(S_d(1:k-1)); % I is the indices where it has the min
            nc(D) = 1;                   
        end
        
        
        % Reset the algorithm
        start = k;
        S_i(start-1) = 0;
        G_i(start-1) = 0; 
        S_d(start-1) = 0;
        G_d(start-1) = 0;
        
    end
    
end
% get where changes are detected
nc_idx=find(nc==1); % estimation change location
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
x_len = length(y); 


false_negative = length(time_detection) - true_positive;
false_negative_rate(iter) = (false_negative)/(true_positive + false_negative);

true_positive_rate(iter)  = true_positive/(true_positive + false_negative);
true_negative = x_len - true_positive;
true_negative_rate(iter) = true_negative/(true_negative + false_positive);
false_positive_rate(iter) = false_positive/(false_positive + true_negative);
end

FPR = mean(false_positive_rate) 
FNR = mean(false_negative_rate)
delayR = mean(delay) 

timeElapsed = toc % to calculat execution time in seconds

%% plot all results
fig = figure;
plot(y),xlim([0,length(y)]),grid on;
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
title('CUSUM Detection in Signal x[n]');
xlabel('n'),ylabel('x[n]');
legend('Location','Southoutside');
legend([h1, h2], {'Actual', 'TP'});

% add title to the figure window
set(gcf,'name','Baseline 2 with an increasing window','numbertitle','off');
savefig(fig,'baseline2_increasing_2d.fig');

% save experiments results and set-up in mat file
save('baseline2_increasing_2d.mat','iterations','channels_num','segments_num','variance','correlation_value',...
'h','d','threshold_detection','x','time_detection','Detection','false_negative',...
'false_positive','true_positive','true_negative','delay','false_positive_rate',...
    'false_negative_rate','delayR','FPR','FNR');
