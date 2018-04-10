%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Yasmin Fathy (Fathy.Yasmin@gmail.com)
% Last Updated: 7 August, 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% The main idea is to generate correlated/uncorrelated multi-channel 
% signals based on a correlation value between them, 
% the no of segments, the no of channels and a variance value. However, we assume that all channels
% will have the same variance
% Test function
% [x, corr_x, time_detection] = randomPieceWiseGenerator(12, 10, 1, 0.5, 1);

function [x, corr_x, time_detection] = randomPieceWiseGenerator(channels_num, segments_num, var, correlation_value,plot_flag)

length_segments = zeros(segments_num, 1); 
segment_len_range1 = 100; 
segment_len_range2 = 500; 
mu_range1 = -3; 
mu_range2 = 3; 
magnitude_range1 = 1; 
magnitude_range2 = 3; 
x = zeros(channels_num, segments_num);
mu = zeros(channels_num, segments_num);

for seg = 1: segments_num
    length_segments(seg) =  randi([segment_len_range1 segment_len_range2],1,1);
    time_detection(seg) = sum(length_segments);
    
    if seg ==1
        mu1 = mu_range1 + (mu_range2 + mu_range2) * rand(1);
        % assign the first segment for all channels the same mu
        mu(:,seg) = repmat(mu1,[seg,1]);
    else
        selected_direction = randi(2);
        for ch = 1: channels_num
            magnitude_value = magnitude_range1 + rand(1,1) * (magnitude_range2-magnitude_range1);
            if selected_direction == 1
                mu(ch,seg) =  mu(ch,seg-1) - magnitude_value;
            else 
                mu(ch,seg) =  mu(ch,seg-1)+ magnitude_value;
            end 
        end
    end
end
start = 1; 
for seg = 1: segments_num
    seg_len = length_segments(seg); 
    endd = start+seg_len-1; 
    
    for ch = 1: channels_num
        x(ch,(start:endd)) = normrnd(mu(ch,seg),var,1,seg_len);
        
     end
    start = endd+1;    
end        

if correlation_value > 0
    correlation_matrix(1:channels_num, 1:channels_num)=correlation_value;
    % assign the diagnoal in the matrix by 1 to create the correlation matrix
    correlation_matrix(logical(eye(size(correlation_matrix)))) = 1;
    C = chol(correlation_matrix);
    corr_x = x' * C;
    corr_x = corr_x';
else
    corr_x = x;
end


time_detection = time_detection(1:end-1);

%% plot and save the result
if plot_flag == 1
    
    fig_g = figure;
    subplot(211),
    % plot uncorrelated signal
    plot(x'),
    ylim = get(gca,'YLim');
    % plot the changes positions on the uncorrelated signal
    for i=1:length(time_detection)
        line([time_detection(i) time_detection(i)], ylim, 'Color', [1 0 0])
    end
    title('Multi-channel Signal x[n]'),
    xlabel('n'),ylabel('x[n]'),
    grid on;
    hold on
    subplot(212),
    % plot the correlated signal
    plot(corr_x');
    ylim = get(gca,'YLim');
    % plot the changes positions on the correlated signal
    for i=1:length(time_detection)
        line([time_detection(i) time_detection(i)], ylim, 'Color', [1 0 0])
    end
    hold off
    grid on;
    title(['Correlated Multi-channel Signal x[n] with corr=' num2str(correlation_value)]);
    xlabel('n'),ylabel('x[n]');
    savefig(fig_g,'randomPieceWiseGenerator.fig');
    %plot(mu);
    
    % save x in a separate figure
    fig = figure;
    % plot uncorrelated signal
    plot(x'),
    ylim = get(gca,'YLim');
    % plot the changes positions on the uncorrelated signal
    for i=1:length(time_detection)
        line([time_detection(i) time_detection(i)], ylim,'Color', [1 0 0])
    end
    xlabel('Samples'),ylabel('Values'),
    grid on;
    savefig(fig,'randomPieceWiseGenerator.fig');    
end
end
