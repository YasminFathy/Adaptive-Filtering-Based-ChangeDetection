

% participant 1
load Participant1 
p1 = Participant1;
[rows, columns] = size(Participant1);

start_p = 9001;
end_p = rows;

x_acc1 = p1{start_p:end_p,1};
y_acc1 = p1{start_p:end_p,2};
z_acc1 = p1{start_p:end_p,3};
labels1= p1{start_p:end_p,4};
vector_magnitude1 = sqrt(x_acc1.^2 + y_acc1.^2 + z_acc1.^2);

[activity_ID1, activity_string] = grp2idx(labels1);
activity_string
[lgt first last e] = SplitVec(activity_ID1, [], 'length','first','last', 'firstelem');

first
last

% walking      1: 9000

% standing    9001: 18000
% jogging    18001: 27000
% sitting    27001 : 36000
% biking     36001 : 45000
% upstairs   45001 : 54000
% downstairs 54001 : 63000

% participant 2
load Participant2
p2 = Participant2;

x_acc2 =  p2{start_p:end_p,1};
y_acc2 =  p2{start_p:end_p,2};
z_acc2 =  p2{start_p:end_p,3};
labels2 = p2{start_p:end_p,4};
vector_magnitude2 = sqrt(x_acc2.^2 + y_acc2.^2 + z_acc2.^2);



% participant 3
load Participant3
p3 = Participant3;

x_acc3 =  p3{start_p:end_p,1};
y_acc3 =  p3{start_p:end_p,2};
z_acc3 =  p3{start_p:end_p,3};
labels3 = p3{start_p:end_p,4};
vector_magnitude3 = sqrt(x_acc3.^2 + y_acc3.^2 + z_acc3.^2);


% participant 4
load Participant4
p4 = Participant4;

x_acc4  = p4{start_p:end_p,1};
y_acc4  = p4{start_p:end_p,2};
z_acc4  = p4{start_p:end_p,3};
labels4 = p4{start_p:end_p,4};
vector_magnitude4 = sqrt(x_acc4.^2 + y_acc4.^2 + z_acc4.^2);



% participant 5
load Participant5
p5 = Participant5;

x_acc5 =  p5{start_p:end_p,1};
y_acc5 =  p5{start_p:end_p,2};
z_acc5 =  p5{start_p:end_p,3};
labels5 = p5{start_p:end_p,4};
vector_magnitude5 = sqrt(x_acc5.^2 + y_acc5.^2 + z_acc5.^2);



plot(vector_magnitude1);
hold on
plot(vector_magnitude2);
plot(vector_magnitude3);
plot(vector_magnitude4);
plot(vector_magnitude5);

ylim = get(gca,'YLim');

 for i=2:length(first)
        %line([time_detection(i) time_detection(i)], ylim, 'LineStyle','--','Color', [1 0 0])
        line([first(i) first(i)], ylim,'Color', [1 0 0])
 end
 
 % combine all semsor data
 % columns will be sensors
x=[vector_magnitude1,vector_magnitude2, vector_magnitude3,... 
    vector_magnitude4, vector_magnitude5];
% transpose to have rows = sensors
x = x';

time_detection = first(2:end);
sensors_num = 5;
