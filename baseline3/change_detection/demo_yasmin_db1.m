%clear;
%%
%%%% Added by Yasmin to have a fair comparion;using the same 
% generated function


iter =1;
iterations = 1;
threshold_detection = 50;

run('real-world-dataset/dataset1.m');

tic % to calculate the execution time in seconds
for iter = 1:iterations
    iter

    y = x;

%%%% Finish updated by Yasmin

%load logwell.mat

% choice of alpha
alpha = 0.01;%.01;

n = 50;  
k = 10; 

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

nd_idx = originalLocations;

% counter
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

%% plots figure

fig = figure;
subplot(2,1,1);

plot(y');
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

% just plot points detected by algorithm
hold on;
ylim = get(gca,'YLim');


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



