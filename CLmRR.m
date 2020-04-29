    % author lfkrieger
% adjusted to experimental group CLm13-18
% start adjustments 04/2020

%% Part1
%% read-in .csv file
% from DeepLabCut output
% coordinates of mouse landmarks (less than with CLm1-12)
% 1-3 nose, 5-7 right ear, 8-10 left ear, 11-13 bodymass, 14-16 redLED_OFF,
% 17-19 blueLED_OFF, 20-22 redLED_ON, 23-25 blueLED_ON
file ='CLm3_20191109_CLm3_20191109_Run1DeepCut_resnet50_CLmicenov28shuffle1_1030000.csv';
T = readtable(file);
T = table2array(T(3:end,:));

%%find direction TS
%red LED likelihood
    Tredon = str2double([T(:,31)]);  %CLm1-10: 31, CLm13-18: 22
    Tredoff = str2double([T(:,25)]); %CLm1-10: 25, CLm13-18: 16

    for i = 1:length(Tredon)-1
        l(i) = (Tredon(i+1,1)-Tredon(i,1));    
    end
    RedOn = find(l>0.95); %frame on %+1 frame
    RedOff = find(l<-0.95); %frame off

%blue LED 
    Tblueon = str2double([T(:,34)]);  %CLm1-10: 34, CLm13-18: 25
    Tblueoff = str2double([T(:,28)]); %CLm1-10: 28, CLm13-18: 19

    for i = 1:length(Tblueon)-1
        l(i) = (Tblueon(i+1,1)-Tblueon(i,1));    
    end
    BlueOn = find(l>0.95); %frame on %+1 frame
    BlueOff = find(l<-0.95); %frame off
    clear l

% bodymass xy plot for Run   
    plot(str2double([T(1:end,11)])), hold on;
    plot(str2double([T(1:end,12)])), title('body mass'), legend('x coord','y coord'), hold on;
    vline(RedOn,'r') %redLED TS
    vline(BlueOn,'b') %blueLED TS

%% manually define TS identity accoring to order in time
% 0 = start/stop
% w/o laser
    % fwbw = 1; bwfw = 2
    % fwfw = -1; bwbw = -2 
% w laser
    % fwbw = 3; bwfw = 4
    % fwfw = -3; bwbw = -4 
% normally sin5: TSid = [0,1,2,1,2,1,2,1,0];
TSid = [0,1,2,1,2,1,2,1,2];

% join red+blue LED TS
    dTS = [RedOn BlueOn];
    dTS = sort(dTS);
    dTS = vertcat(dTS,TSid);
    dTS = dTS';
    
%%
T = readtable(file); %because I want the original table be saved 
m18.D20200301.Run5={T,dTS};
clear ans BlueOn BlueOff dTS i RedOn RedOff T Tblueon Tblueoff Tredon Tredoff TSid file

% m15.D20200227.Run5{1}=T;  

%% Part2 
%% post-processing
%explanation of data structure:
%strct containing table
% to make calculation, convert to double: table2array (tabel to cell),
% str2double (cell to double)
% once done, convert back to original strucutre: num2cell(double to cell), then cell2table (cell to table) 

%1.) normalize coordinates to zero position = mouse sitting still before
%run onset or in first 40*9 frames (9secs @40Hz)
% = transform to relative mouse displacement in pixels from resting-position

fns = sort(fieldnames(m7));
for i = 1:length(fns) % for all days
    fns2 = fieldnames(m7.(fns{i}));
    vec = [2,3,5,6,8,9,11,12,14,15,17,18,20,21,23,24,26,27,29,30,32,33]; %for CLm1-10
    %vec = [2,3,5,6,8,9,11,12,14,15,17,18,20,21,23,24]; %for CLm13-18
    for j = 1:length(fns2) % for all Runs
        T = table2array(m7.(fns{i}).(fns2{j}){1,1}(1:end,:));
        for k = 1:length(vec)% for all x y columns
            %normalize to mean of first 150 frames (CLm1-10)/360 frames (CLm13-18)
            M = mean(str2double(table2array(m7.(fns{i}).(fns2{j}){1,1}(3:152,vec(k))))); %mean of first 360 frames
            T(3:end,vec(k)) = num2cell(str2double(table2array(m7.(fns{i}).(fns2{j}){1,1}(3:end,vec(k))))-M);  %normalize
        end
        m7.(fns{i}).(fns2{j}){1,1}=cell2table(T);
    end
end

%2.) excluding coord with low estimation-likelihood =NaN
% but what cutoff value? 95%
fns = fieldnames(m7);
for i = 1:length(fns) % for all days
    fns2 = fieldnames(m7.(fns{i}));
    for j = 1:length(fns2) % for all Runs
        T = (table2array(m7.(fns{i}).(fns2{j}){1,1}(1:end,:)));
        idx=find(str2double(T(:,7))<0.95);
        for k=1:length(idx)
            T{idx(k),5}=NaN; 
            T{idx(k),6}=NaN; %for right ear
        end
        idx=find(str2double(T(:,19))<0.95);
        for k=1:length(idx)
            T{idx(k),17}=NaN;
            T{idx(k),18}=NaN; %for center dody mass
        end
        %dummy2=(str2double(table2array(m7.D20191107.Run4{1,1}(dummy,7))));
        %histogram(dummy2,50)
        m7.(fns{i}).(fns2{j}){1,1}=cell2table(T);
    end
end
clearvars -except m7 

%% get all TS of a type and associated timeseries of landmark coord.

%here for two landmarks, the center of bodymass and the right ear
l_fwbw =cell(4,2); %{1,1} = bodymass x, {1,2} = bodymass y
l_bwfw =cell(4,2); %{1,1} = bodymass x, {1,2} = bodymass y
e_fwbw =cell(4,2); %{1,1} = right ear x, {1,2} = right ear y
e_bwfw =cell(4,2); %{1,1} = right ear x, {1,2} = right ear y

fns = fieldnames(m7);
for i = 1:length(fns) % for all days
    fns2 = fieldnames(m7.(fns{i}));
    for j = 1:length(fns2) % for all Runs
    % find all TS of certain kind + associated frame number    
        %fwbw w/o laser
        idx_fw = find(m7.(fns{i}).(fns2{j}){1,2}(:,2) ==1);
        idx_fw2 = m7.(fns{i}).(fns2{j}){1,2}(idx_fw,1)+3;
        %bwfw w/o laser
        idx_bw = find(m7.(fns{i}).(fns2{j}){1,2}(:,2) ==2);
        idx_bw2 = m7.(fns{i}).(fns2{j}){1,2}(idx_bw,1)+3; 
        %fwbw w/ laser
        idx_fwL = find(m7.(fns{i}).(fns2{j}){1,2}(:,2) ==3);
        idx_fwL2 = m7.(fns{i}).(fns2{j}){1,2}(idx_fwL,1)+3;
        %bwfw w/laser
        idx_bwL = find(m7.(fns{i}).(fns2{j}){1,2}(:,2) ==4);
        idx_bwL2 = m7.(fns{i}).(fns2{j}){1,2}(idx_bwL,1)+3;
        %fwfw w/o laser
        idx_fwfw = find(m7.(fns{i}).(fns2{j}){1,2}(:,2) ==-1);
        idx_fwfw2 = m7.(fns{i}).(fns2{j}){1,2}(idx_fwfw,1)+3;
        %bwbw w/o laser
        idx_bwbw = find(m7.(fns{i}).(fns2{j}){1,2}(:,2) ==-2);
        idx_bwbw2 = m7.(fns{i}).(fns2{j}){1,2}(idx_bwbw,1)+3;
        %fwfw w/ laser
        idx_fwfwL = find(m7.(fns{i}).(fns2{j}){1,2}(:,2) ==-3);
        idx_fwfwL2 = m7.(fns{i}).(fns2{j}){1,2}(idx_fwfwL,1)+3;
        %bwbw w/ laser
        idx_bwbwL = find(m7.(fns{i}).(fns2{j}){1,2}(:,2) ==-4);
        idx_bwbwL2 = m7.(fns{i}).(fns2{j}){1,2}(idx_bwbwL,1)+3;
    % find tracking coordinates of body mass and right ear around TS
        T = table2array(m7.(fns{i}).(fns2{j}){1,1}(1:end,:));
        %l_fwbw AND e_fwbw
        for k = 1:length(idx_fw2) %normal fwbw x and y coord 
            l_fwbw{1,1}(:,size(l_fwbw{1,1},2)+1)= (T(idx_fw2(k)-40*6:idx_fw2(k)+40*6,17)); % bodymass x  %17 for CLm1-10, 11 for CLm13-18
            l_fwbw{1,2}(:,size(l_fwbw{1,2},2)+1)= (T(idx_fw2(k)-40*6:idx_fw2(k)+40*6,18)); % bodymass y  %18 for CLm1-10, 12 for CLm13-18
            e_fwbw{1,1}(:,size(e_fwbw{1,1},2)+1)= (T(idx_fw2(k)-40*6:idx_fw2(k)+40*6,5)); % right ear x
            e_fwbw{1,2}(:,size(e_fwbw{1,2},2)+1)= (T(idx_fw2(k)-40*6:idx_fw2(k)+40*6,6)); % right ear y
        end
        clear k
        for k = 1:length(idx_fwL2) %normal fwbw x and y coord WITH LASER 
            l_fwbw{2,1}(:,size(l_fwbw{2,1},2)+1)= (T(idx_fwL2(k)-40*6:idx_fwL2(k)+40*6,17)); % bodymass x
            l_fwbw{2,2}(:,size(l_fwbw{2,2},2)+1)= (T(idx_fwL2(k)-40*6:idx_fwL2(k)+40*6,18)); % bodymass y
            e_fwbw{2,1}(:,size(e_fwbw{2,1},2)+1)= (T(idx_fwL2(k)-40*6:idx_fwL2(k)+40*6,5)); % right ear x
            e_fwbw{2,2}(:,size(e_fwbw{2,2},2)+1)= (T(idx_fwL2(k)-40*6:idx_fwL2(k)+40*6,6)); % right ear y
        end
        clear k
        for k = 1:length(idx_fwfw2) %flex fwfw x and y coord 
            l_fwbw{3,1}(:,size(l_fwbw{3,1},2)+1)= (T(idx_fwfw2(k)-40*6:idx_fwfw2(k)+40*6,17)); % bodymass x
            l_fwbw{3,2}(:,size(l_fwbw{3,2},2)+1)= (T(idx_fwfw2(k)-40*6:idx_fwfw2(k)+40*6,18)); % bodymass y
            e_fwbw{3,1}(:,size(e_fwbw{3,1},2)+1)= (T(idx_fwfw2(k)-40*6:idx_fwfw2(k)+40*6,5)); % right ear x
            e_fwbw{3,2}(:,size(e_fwbw{3,2},2)+1)= (T(idx_fwfw2(k)-40*6:idx_fwfw2(k)+40*6,6)); % right ear y
        end
        clear k
        for k = 1:length(idx_fwfwL2) %flex fwfw x and y coord WITH LASER     
            l_fwbw{4,1}(:,size(l_fwbw{4,1},2)+1)= (T(idx_fwfwL2(k)-40*6:idx_fwfwL2(k)+40*6,17)); % bodymass x
            l_fwbw{4,2}(:,size(l_fwbw{4,2},2)+1)= (T(idx_fwfwL2(k)-40*6:idx_fwfwL2(k)+40*6,18)); % bodymass y
            e_fwbw{4,1}(:,size(e_fwbw{4,1},2)+1)= (T(idx_fwfwL2(k)-40*6:idx_fwfwL2(k)+40*6,5)); % right ear x
            e_fwbw{4,2}(:,size(e_fwbw{4,2},2)+1)= (T(idx_fwfwL2(k)-40*6:idx_fwfwL2(k)+40*6,6)); % right ear y
        end
        clear k    
        % l_bwfw
        for k = 1:length(idx_bw2) %normal bwfw x and y coord
            l_bwfw{1,1}(:,size(l_bwfw{1,1},2)+1)= (T(idx_bw2(k)-40*6:idx_bw2(k)+40*6,17)); % bodymass x
            l_bwfw{1,2}(:,size(l_bwfw{1,2},2)+1)= (T(idx_bw2(k)-40*6:idx_bw2(k)+40*6,18)); % bodymass y
            e_bwfw{1,1}(:,size(e_bwfw{1,1},2)+1)= (T(idx_bw2(k)-40*6:idx_bw2(k)+40*6,5)); % right ear x
            e_bwfw{1,2}(:,size(e_bwfw{1,2},2)+1)= (T(idx_bw2(k)-40*6:idx_bw2(k)+40*6,6)); % right ear y
        end
        clear k
        for k = 1:length(idx_bwL2) %normal bwfw x and y coord WITH LASER
            l_bwfw{2,1}(:,size(l_bwfw{2,1},2)+1)= (T(idx_bwL2(k)-40*6:idx_bwL2(k)+40*6,17)); % bodymass x
            l_bwfw{2,2}(:,size(l_bwfw{2,2},2)+1)= (T(idx_bwL2(k)-40*6:idx_bwL2(k)+40*6,18)); % bodymass y
            e_bwfw{2,1}(:,size(e_bwfw{2,1},2)+1)= (T(idx_bwL2(k)-40*6:idx_bwL2(k)+40*6,5)); % right ear x
            e_bwfw{2,2}(:,size(e_bwfw{2,2},2)+1)= (T(idx_bwL2(k)-40*6:idx_bwL2(k)+40*6,6)); % right ear y
        end
        clear k
        for k = 1:length(idx_bwbw2) %flex bwbw x and y coord
            l_bwfw{3,1}(:,size(l_bwfw{3,1},2)+1)= (T(idx_bwbw2(k)-40*6:idx_bwbw2(k)+40*6,17)); % bodymass x
            l_bwfw{3,2}(:,size(l_bwfw{3,2},2)+1)= (T(idx_bwbw2(k)-40*6:idx_bwbw2(k)+40*6,18)); % bodymass y
            e_bwfw{3,1}(:,size(e_bwfw{3,1},2)+1)= (T(idx_bwbw2(k)-40*6:idx_bwbw2(k)+40*6,5)); % rigt ear x
            e_bwfw{3,2}(:,size(e_bwfw{3,2},2)+1)= (T(idx_bwbw2(k)-40*6:idx_bwbw2(k)+40*6,6)); % right ear y
        end
        clear k
        for k = 1:length(idx_bwbwL2) %flex bwbw x and y coord WITH LASER
            l_bwfw{4,1}(:,size(l_bwfw{4,1},2)+1)= (T(idx_bwbwL2(k)-40*6:idx_bwbwL2(k)+40*6,17)); % bodymass x
            l_bwfw{4,2}(:,size(l_bwfw{4,2},2)+1)= (T(idx_bwbwL2(k)-40*6:idx_bwbwL2(k)+40*6,18)); % bodymass y
            e_bwfw{4,1}(:,size(e_bwfw{4,1},2)+1)= (T(idx_bwbwL2(k)-40*6:idx_bwbwL2(k)+40*6,5)); % right ear x
            e_bwfw{4,2}(:,size(e_bwfw{4,2},2)+1)= (T(idx_bwbwL2(k)-40*6:idx_bwbwL2(k)+40*6,6)); % right ear y
        end
        clear k             
    end  
end
clearvars -except m7 l_bwfw l_fwbw e_bwfw e_fwbw

%% plot landmark coordinates aligned to direction TS

% here has to be done for each landmark individually
% raw values
    var = l_fwbw;
    var2 = 'fwbw bodymass x';
    var3 = 'fwbw bodymass y';
        figure %fwbw
        subplot(6,2,1)%fwbw RR protocol
        x=0:480; y=linspace(16,18,80); y=[y linspace(5,5,160)]; y=[y linspace(-10,-16,241)]; y2=NaN(1,240); y2=[y2 linspace(10,16,241)]; 
        plot(x,y,'color', [0.75 0.75 0.75]), hold on, plot(x,y2,'r'), vline([80 241],{'b','b'}), hline(0,'k'), xlim([0 481]), set(gca,'XTick',[],'YTick',[-15 15]); ylabel('rpm');
        subplot(6,2,2)
        plot(x,y,'color', [0.75 0.75 0.75]), hold on, plot(x,y2,'r'),vline([80 241],{'b','b'}), hline(0,'k'), xlim([0 481]), set(gca,'XTick',[],'YTick',[-15 15]); ylabel('rpm');
        subplot(6,2,[3 5 7 9 11]) %x coord
        p1=plot(cell2mat(var{1,1}),':', 'color', [0.75 0.75 0.75]); hold on
        p2=plot(cell2mat(var{2,1}),'b','LineWidth',1.5); hold on
        p3=plot(cell2mat(var{3,1}), 'g','LineWidth',1.5); hold on
        p4=plot(cell2mat(var{4,1}), 'r','LineWidth',1.5);
        p5=plot(nanmean(cell2mat(var{1,1})'),'color', [0.2 0.2 0.2],'LineWidth',2);
        xlim([0 481]), ylim([-40 40]), vline([80 241],{'b','b'}), title(var2), ylabel('lateral coord.'), xlabel('[seconds]');
        ax = gca;    
        ax.XTick = [0 40 80 120 160 200 240 280 320 360 400 440 480];        
        ax.XTickLabel = {'-6','-5','-4','-3','-2','-1','0','1','2','3','4','5','6'};    
        subplot(6,2,[4 6 8 10 12]) % y coord
        plot(cell2mat(var{1,2}),':', 'color', [0.75 0.75 0.75]), hold on
        plot(cell2mat(var{2,2}), 'b','LineWidth',1.5), hold on
        plot(cell2mat(var{3,2}), 'g','LineWidth',1.5), hold on
        plot(cell2mat(var{4,2}), 'r','LineWidth',1.5),
        plot(nanmean(cell2mat(var{1,2})'),'color', [0.2 0.2 0.2],'LineWidth',2);
        xlim([0 481]), ylim([-100 80]), vline([80 241],{'b','b'}), title(var3), ylabel('longitudinal coord.'), xlabel('[seconds]');
        ax = gca;    
        ax.XTick = [0 40 80 120 160 200 240 280 320 360 400 440 480];        
        ax.XTickLabel = {'-6','-5','-4','-3','-2','-1','0','1','2','3','4','5','6'};
        legend([p1(1) p2(1) p3(1) p4(1) p5(1)],'normal','normal+laser','flex','flex+laser','mean normal')  
     
    var = l_bwfw;
    var2 = 'bwfw bodymass x';
    var3 = 'bwfw bodymass y';
        figure %bwfw
        subplot(6,2,1)%bwfw RR protocol
        x=0:480; y=linspace(-16,-18,80); y=[y linspace(-5,-5,160)]; y=[y linspace(10,16,241)]; y2=NaN(1,240); y2=[y2 linspace(-10,-16,241)]; 
        plot(x,y), hold on, plot(x,y2,'r'), vline([80 241],{'b','b'}), hline(0,'k'), xlim([0 481]), set(gca,'XTick',[],'YTick',[-15 15]); ylabel('rpm');
        subplot(6,2,2)
        plot(x,y), hold on, plot(x,y2,'r'), vline([80 241],{'b','b'}), hline(0,'k'), xlim([0 481]), set(gca,'XTick',[],'YTick',[-15 15]); ylabel('rpm');
        subplot(6,2,[3 5 7 9 11]) %x coord
        p1=plot(cell2mat(var{1,1}),':', 'color', [0.75 0.75 0.75]); hold on
        p2=plot(cell2mat(var{2,1}), 'b','LineWidth',1.5); hold on
        p3=plot(cell2mat(var{3,1}), 'g','LineWidth',1.5); hold on
        p4=plot(cell2mat(var{4,1}), 'r','LineWidth',1.5);
        p5=plot(nanmean(cell2mat(var{1,1})'),'color', [0.2 0.2 0.2],'LineWidth',2);
        xlim([0 481]),ylim([-40 40]), vline([80 241],{'b','b'}), title(var2), ylabel('lateral coord.'), xlabel('[seconds]');
        ax = gca;    
        ax.XTick = [0 40 80 120 160 200 240 280 320 360 400 440 480];        
        ax.XTickLabel = {'-6','-5','-4','-3','-2','-1','0','1','2','3','4','5','6'};
        subplot(6,2,[4 6 8 10 12]) %y coord
        plot(cell2mat(var{1,2}),':', 'color', [0.75 0.75 0.75]), hold on
        plot(cell2mat(var{2,2}), 'b','LineWidth',1.5), hold on
        plot(cell2mat(var{3,2}), 'g','LineWidth',1.5), hold on
        plot(cell2mat(var{4,2}), 'r','LineWidth',1.5)
        plot(nanmean(cell2mat(var{1,2})'),'color', [0.2 0.2 0.2],'LineWidth',2);
        xlim([0 481]), ylim([-50 60]), vline([80 241],{'b','b'}), title(var3), ylabel('longitudinal coord.'), xlabel('[seconds]');
        ax = gca;    
        ax.XTick = [0 40 80 120 160 200 240 280 320 360 400 440 480];        
        ax.XTickLabel = {'-6','-5','-4','-3','-2','-1','0','1','2','3','4','5','6'};
        legend([p1(1) p2(1) p3(1) p4(1) p5(1)],'normal','normal+laser','flex','flex+laser','mean normal')  
       
%% z-transformation and plot

% ztransformation of the whole cell array _fwbw or _bwfw
 var = l_bwfw; %! RENAME var_z output variable!   
    [rows cols] = size(var);
    for g = 1:cols
        % mean and sd of normal TS {1,1} and {1,2} row-wise (of all reps)
        for gg= 1:length(var{1,g}(:,1))
            m(gg,1) = nanmean(cell2mat(var{1,g}(gg,:)));
            s(gg,1) =  nanstd(cell2mat(var{1,g}(gg,:)));
        end
        for h = 1:rows
            if isempty(cell2mat(var{h,g})) %here cell2mat?
                continue
            end
            for j=1:length(var{h,g}(:,1)) %for rows
                for i = 1:length(var{h,g}(1,:)) %for each column value (rep)
                   var_z{h,g}(j,i)=(cell2mat(var{h,g}(j,i))-m(j))/s(j);
                end
            end
        end
    end
    clearvars g gg h j i cols rows s m yl gray ans var

% plot ztransformed coord of repetitions per timepoint 
% !must change VAR: landmark and fwbw!
    var = l_fwbw_z;
    var2 = 'fwbw bodymass x';
    var3 = 'fwbw bodymass y';
        figure %fwbw
        subplot(6,2,1)%fwbe RR protocol
        x=0:480; y=linspace(16,18,80); y=[y linspace(5,5,160)]; y=[y linspace(-10,-16,241)]; y2=NaN(1,240); y2=[y2 linspace(10,16,241)]; 
        plot(x,y,'color', [0.75 0.75 0.75]), hold on, plot(x,y2,'r'), vline([80 241],{'b','b'}), hline(0,'k'), xlim([0 481]), set(gca,'XTick',[],'YTick',[-15 15]); ylabel('rpm');
        subplot(6,2,2)
        plot(x,y,'color', [0.75 0.75 0.75]), hold on, plot(x,y2,'r'),vline([80 241],{'b','b'}), hline(0,'k'), xlim([0 481]), set(gca,'XTick',[],'YTick',[-15 15]); ylabel('rpm');
        subplot(6,2,[3 5 7 9 11]) %x coord
        [rows cols] = size(var{1,1});
        p1=plot((var{1,1}),':', 'color', [0.75 0.75 0.75]); hold on
        p2=plot((var{2,1}), 'b','LineWidth',1.5); hold on
        p3=plot((var{3,1}), 'g','LineWidth',1.5); hold on
        p4=plot((var{4,1}), 'r','LineWidth',1.5);
        line([0 rows],[2 2],'Color',[0 0 0])
        line([0 rows],[-2 -2],'Color',[0 0 0])
        xlim([0 481]), ylim([-5 5]), vline([80 241],{'b','b'}), title(var2), ylabel('ztransformed lateral coord.'), xlabel('[seconds]');
        ax = gca;    
        ax.XTick = [0 40 80 120 160 200 240 280 320 360 400 440 480];        
        ax.XTickLabel = {'-6','-5','-4','-3','-2','-1','0','1','2','3','4','5','6'};
        subplot(6,2,[4 6 8 10 12]) % y coord
        [rows cols] = size(var{1,2});
        plot((var{1,2}),':', 'color', [0.75 0.75 0.75]); hold on
        plot((var{2,2}), 'b','LineWidth',1.5); hold on
        plot((var{3,2}), 'g','LineWidth',1.5); hold on
        plot((var{4,2}), 'r','LineWidth',1.5);
        line([0 rows],[2 2],'Color',[0 0 0])
        line([0 rows],[-2 -2],'Color',[0 0 0])
        xlim([0 481]), ylim([-7 5]), vline([80 241],{'b','b'}), title(var3), ylabel('ztransformed longitudinal coord.'), xlabel('[seconds]');
        ax = gca;    
        ax.XTick = [0 40 80 120 160 200 240 280 320 360 400 440 480];        
        ax.XTickLabel = {'-6','-5','-4','-3','-2','-1','0','1','2','3','4','5','6'};
        legend([p1(1) p2(1) p3(1) p4(1)],'normal','normal+laser','flex', 'flex+laser')  
 
    var = e_bwfw_z;
    var2 = 'bwfw right ear x';
    var3 = 'bwfw right ear y';
        figure %bwfw
        subplot(6,2,1)%bwfw RR protocol
        x=0:480; y=linspace(-16,-18,80); y=[y linspace(-5,-5,160)]; y=[y linspace(10,16,241)]; y2=NaN(1,240); y2=[y2 linspace(-10,-16,241)]; 
        plot(x,y), hold on, plot(x,y2,'r'), vline([80 241],{'b','b'}), hline(0,'k'), xlim([0 481]), set(gca,'XTick',[],'YTick',[-15 15]); ylabel('rpm');
        subplot(6,2,2)
        plot(x,y), hold on, plot(x,y2,'r'), vline([80 241],{'b','b'}), hline(0,'k'), xlim([0 481]), set(gca,'XTick',[],'YTick',[-15 15]); ylabel('rpm');
        subplot(6,2,[3 5 7 9 11]) %x coord
        [rows cols] = size(var{1,1});
        p1=plot((var{1,1}),':', 'color', [0.75 0.75 0.75]); hold on
        p2=plot((var{2,1}), 'b','LineWidth',1.5); hold on
        p3=plot((var{3,1}), 'g','LineWidth',1.5); hold on
        p4=plot((var{4,1}), 'r','LineWidth',1.5);
        line([0 rows],[2 2],'Color',[0 0 0])
        line([0 rows],[-2 -2],'Color',[0 0 0])
        xlim([0 481]), ylim([-5 5]), vline([80 241],{'b','b'}), title(var2), ylabel('ztransformed lateral coord.'), xlabel('[seconds]');
        ax = gca;    
        ax.XTick = [0 40 80 120 160 200 240 280 320 360 400 440 480];        
        ax.XTickLabel = {'-6','-5','-4','-3','-2','-1','0','1','2','3','4','5','6'};
        subplot(6,2,[4 6 8 10 12]) % y coord
        [rows cols] = size(var{1,2});
        plot((var{1,2}),':', 'color', [0.75 0.75 0.75]), hold on
        plot((var{2,2}), 'b','LineWidth',1.5), hold on
        plot((var{3,2}), 'g','LineWidth',1.5), hold on
        plot((var{4,2}), 'r','LineWidth',1.5)
        line([0 rows],[2 2],'Color',[0 0 0])
        line([0 rows],[-2 -2],'Color',[0 0 0])
        xlim([0 481]), ylim([-5 5]), vline([80 241],{'b','b'}), title(var3), ylabel('ztransformed longitudinal coord.'), xlabel('[seconds]'); 
        ax = gca;    
        ax.XTick = [0 40 80 120 160 200 240 280 320 360 400 440 480];        
        ax.XTickLabel = {'-6','-5','-4','-3','-2','-1','0','1','2','3','4','5','6'};
        legend([p1(1) p2(1) p3(1) p4(1)],'normal','normal+laser','flex','flex+laser')  
        
%% All mice analysis
%1.) create a table containing all TS-aligned snippets

    MouseID = {'m1','m3','m6','m7','m14','m15','m16','m17','m18'}';
    ArchT   = logical([1;1;0;0;0;1;1;0;1]);
    Control = logical([0;0;1;1;1;0;0;1;0]);
    %mX_file = {m1 m3 m6 m7 m14 m15 m16 m17 m18}';
    l_fwbw  = {m1_l_fwbw m3_l_fwbw m6_l_fwbw m7_l_fwbw m14_l_fwbw m15_l_fwbw m16_l_fwbw m17_l_fwbw m18_l_fwbw}';
    l_fwbw_z= {m1_l_fwbw_z m3_l_fwbw_z m6_l_fwbw_z m7_l_fwbw_z m14_l_fwbw_z m15_l_fwbw_z m16_l_fwbw_z m17_l_fwbw_z m18_l_fwbw_z}';
    e_fwbw  = {m1_e_fwbw m3_e_fwbw m6_e_fwbw m7_e_fwbw m14_e_fwbw m15_e_fwbw m16_e_fwbw m17_e_fwbw m18_e_fwbw}';
    e_fwbw_z= {m1_e_fwbw_z m3_e_fwbw_z m6_e_fwbw_z m7_e_fwbw_z m14_e_fwbw_z m15_e_fwbw_z m16_e_fwbw_z m17_e_fwbw_z m18_e_fwbw_z}';
    l_bwfw  = {m1_l_bwfw m3_l_bwfw m6_l_bwfw m7_l_bwfw m14_l_bwfw m15_l_bwfw m16_l_bwfw m17_l_bwfw m18_l_bwfw }';
    l_bwfw_z= {m1_l_bwfw_z m3_l_bwfw_z m6_l_bwfw_z m7_l_bwfw_z m14_l_bwfw_z m15_l_bwfw_z m16_l_bwfw_z m17_l_bwfw_z m18_l_bwfw_z }';
    e_bwfw  = {m1_e_bwfw m3_e_bwfw m6_e_bwfw m7_e_bwfw m14_e_bwfw m15_e_bwfw m16_e_bwfw m17_e_bwfw m18_e_bwfw }';
    e_bwfw_z= {m1_e_bwfw_z m3_e_bwfw_z m6_e_bwfw_z m7_e_bwfw_z m14_e_bwfw_z m15_e_bwfw_z m16_e_bwfw_z m17_e_bwfw_z m18_e_bwfw_z }';
    Tclm = table(MouseID, ArchT, Control, l_fwbw, l_fwbw_z, e_fwbw, e_fwbw_z, l_bwfw, l_bwfw_z, e_bwfw, e_bwfw_z);

% 2.) get mean of means of all mice ArchT normal runs fwbw and bwfw
    for i = 1:9
       if Tclm.ArchT(i) == 1
           m(i,:) = nanmean(cell2mat(Tclm.e_fwbw{i,1}{1,2}),2)'; %row-wise mean       
       else
          continue %how to skip?
       end
    end
    m_e_fwbw=m';
    mm_e_fwbw = mean(m_e_fwbw,2);
    clear m

    for i = 1:9
       if Tclm.ArchT(i) == 1
           m(i,:) = nanmean(cell2mat(Tclm.e_bwfw{i,1}{1,2}),2)'; %row-wise mean       
       else
          continue %how to skip?
       end
    end
    m_e_bwfw=m';
    mm_e_bwfw = mean(m_e_bwfw,2);
    clear m

        plot(mm_e_fwbw, 'k','LineWidth',2), hold on
        plot(mm_e_bwfw,':', 'color','k','LineWidth',2), hold on

% 2.1) get mean of means of all mice ArchT normal runs+laser runs fwbw and bwfw
    for i = 1:9
       if Tclm.ArchT(i) == 1
           m(i,:) = nanmean(cell2mat(Tclm.e_fwbw{i,1}{2,2}),2)'; %row-wise mean       
       else
          continue %how to skip?
       end
    end
    m_b_e_fwbw=m';
    mm_b_e_fwbw = mean(m_b_e_fwbw,2);
    clear m

    for i = 1:9
       if Tclm.ArchT(i) == 1
           m(i,:) = nanmean(cell2mat(Tclm.e_bwfw{i,1}{2,2}),2)'; %row-wise mean       
       else
          continue %how to skip?
       end
    end
    m_b_e_bwfw=m';
    mm_b_e_bwfw = mean(m_b_e_bwfw,2);
    clear m

        plot(mm_b_e_fwbw, 'b','LineWidth',2), hold on
        plot(mm_b_e_bwfw,':', 'color','b','LineWidth',2), hold on

% 2.2) get mean of means of all mice ArchT flex runs fwbw and bwfw
    x=[2 6 7 9];
    for i = 1:4
        j = x(i);
       %if Tclm.ArchT(i) == 1
           m(i,:) = nanmean(cell2mat(Tclm.e_fwbw{j,1}{3,2}),2)'; %row-wise mean       
      % else
      %    continue %how to skip?
      % end
    end
    m_g_e_fwbw=m';
    mm_g_e_fwbw = mean(m_g_e_fwbw,2);
    clear m

    for i = 1:4
        j = x(i);
       %if Tclm.ArchT(i) == 1
           m(i,:) = nanmean(cell2mat(Tclm.e_bwfw{j,1}{3,2}),2)'; %row-wise mean       
       %else
       %   continue %how to skip?
       %end
    end
    m_g_e_bwfw=m';
    mm_g_e_bwfw = mean(m_g_e_bwfw,2);
    clear m

        plot(mm_g_e_fwbw, 'g','LineWidth',2), hold on
        plot(mm_g_e_bwfw,':', 'color','g','LineWidth',2), hold on
        
% 2.3) get mean of means of all mice ArchT flex runs+laser runs fwbw and bwfw
    for i = 1:9
       if Tclm.ArchT(i) == 1
           m(i,:) = nanmean(cell2mat(Tclm.e_fwbw{i,1}{4,2}),2)'; %row-wise mean       
       else
          continue %how to skip?
       end
    end
    m_r_e_fwbw=m';
    mm_r_e_fwbw = mean(m_r_e_fwbw,2);
    clear m

    for i = 1:9
       if Tclm.ArchT(i) == 1
           m(i,:) = nanmean(cell2mat(Tclm.e_bwfw{i,1}{4,2}),2)'; %row-wise mean       
       else
          continue %how to skip?
       end
    end
    m_r_e_bwfw=m';
    mm_r_e_bwfw = mean(m_r_e_bwfw,2);
    clear m

        plot(mm_r_e_fwbw, 'r','LineWidth',2), hold on
        plot(mm_r_e_bwfw,':', 'color','r','LineWidth',2), hold on
        vline([80 241],{'b','b'}), xlim([0 481])
        title('ArchT'), ylabel('displacement along running direction [pxl]'), xlabel('[seconds]'); 
        ax = gca;    
        ax.XTick = [0 40 80 120 160 200 240 280 320 360 400 440 480];        
        ax.XTickLabel = {'-6','-5','-4','-3','-2','-1','0','1','2','3','4','5','6'};

legend('normal','normal','normal+laser','normal+laser')
legend('normal','normal','flex','flex')
legend('normal','normal','flex+laser','flex+laser')
legend('normal','normal','normal+laser','normal+laser','flex','flex','flex+laser','flex+laser')       
     

%% 3.) get mean of means of all mice Control normal runs fwbw and bwfw
    for i = 1:9
       if Tclm.ArchT(i) == 0
           n(i,:) = nanmean(cell2mat(Tclm.e_fwbw{i,1}{1,2}),2)'; %row-wise mean       
       else
          continue %how to skip?
       end
    end
    n_e_fwbw=n';
    nn_e_fwbw = mean(n_e_fwbw,2);
    clear n

    for i = 1:9
       if Tclm.ArchT(i) == 0
           n(i,:) = nanmean(cell2mat(Tclm.e_bwfw{i,1}{1,2}),2)'; %row-wise mean       
       else
          continue %how to skip?
       end
    end
    n_e_bwfw=n';
    nn_e_bwfw = mean(n_e_bwfw,2);
    clear n

        plot(nn_e_fwbw, 'k','LineWidth',2), hold on
        plot(nn_e_bwfw,':', 'color','k','LineWidth',2), hold on

        
% 3.1) get mean of means of all mice Control normal runs+laser runs fwbw and bwfw
    for i = 1:9
       if Tclm.ArchT(i) == 0
           n(i,:) = nanmean(cell2mat(Tclm.e_fwbw{i,1}{2,2}),2)'; %row-wise mean       
       else
          continue %how to skip?
       end
    end
    n_b_e_fwbw=n';
    nn_b_e_fwbw = mean(n_b_e_fwbw,2);
    clear n

    for i = 1:9
       if Tclm.ArchT(i) == 0
           n(i,:) = nanmean(cell2mat(Tclm.e_bwfw{i,1}{2,2}),2)'; %row-wise mean       
       else
          continue %how to skip?
       end
    end
    n_b_e_bwfw=n';
    nn_b_e_bwfw = mean(n_b_e_bwfw,2);
    clear n

        plot(nn_b_e_fwbw, 'b','LineWidth',2), hold on
        plot(nn_b_e_bwfw,':', 'color','b','LineWidth',2), hold on

% 3.2) get mean of means of all mice Control flex runs fwbw and bwfw
    x=[4 5 8];
    for i = 1:3
        j = x(i);
       %if Tclm.ArchT(i) == 1
           n(i,:) = nanmean(cell2mat(Tclm.e_fwbw{j,1}{3,2}),2)'; %row-wise mean       
      % else
      %    continue %how to skip?
      % end
    end
    n_g_e_fwbw=n';
    nn_g_e_fwbw = mean(n_g_e_fwbw,2);
    clear n

    for i = 1:3
        j = x(i);
       %if Tclm.ArchT(i) == 1
           n(i,:) = nanmean(cell2mat(Tclm.e_bwfw{j,1}{3,2}),2)'; %row-wise mean       
       %else
       %   continue %how to skip?
       %end
    end
    n_g_e_bwfw=n';
    nn_g_e_bwfw = mean(n_g_e_bwfw,2);
    clear n

        plot(nn_g_e_fwbw, 'g','LineWidth',2), hold on
        plot(nn_g_e_bwfw,':', 'color','g','LineWidth',2), hold on
        
% 2.3) get mean of means of all mice Control flex runs+laser runs fwbw and bwfw
    for i = 1:9
       if Tclm.ArchT(i) == 0
           n(i,:) = nanmean(cell2mat(Tclm.e_fwbw{i,1}{4,2}),2)'; %row-wise mean       
       else
          continue %how to skip?
       end
    end
    n_r_e_fwbw=n';
    nn_r_e_fwbw = mean(n_r_e_fwbw,2);
    clear n

    for i = 1:9
       if Tclm.ArchT(i) == 0
           n(i,:) = nanmean(cell2mat(Tclm.e_bwfw{i,1}{4,2}),2)'; %row-wise mean       
       else
          continue %how to skip?
       end
    end
    n_r_e_bwfw=n';
    nn_r_e_bwfw = mean(n_r_e_bwfw,2);
    clear n

        plot(nn_r_e_fwbw, 'r','LineWidth',2), hold on
        plot(nn_r_e_bwfw,':', 'color','r','LineWidth',2), hold on
        
        vline([80 241],{'b','b'}), xlim([0 481])
        title('Control'), ylabel('displacement along running direction [pxl]'), xlabel('[seconds]'); 
        ax = gca;    
        ax.XTick = [0 40 80 120 160 200 240 280 320 360 400 440 480];        
        ax.XTickLabel = {'-6','-5','-4','-3','-2','-1','0','1','2','3','4','5','6'};

        
legend('normal fwbw','normal bwfw','normal fwbw+laser','normal bwfw+laser')
legend('normal fwbw ','normal bwfw','flex fwfw','flex bwbw')
legend('normal fwbw','normal bwfw','flex fwfw+laser','flex bwbw+laser')
legend('normal fwbw','normal bwfw','normal fwbw+laser','normal bwfw+laser','flex fwfw','flex bwbw','flex fwfw+laser','flex bwbw+laser')       
          

%%
m18_l_fwbw  = l_fwbw;
m18_l_fwbw_z= l_fwbw_z;
m18_e_fwbw  = e_fwbw;
m18_e_fwbw_z= e_fwbw_z;
m18_l_bwfw  = l_bwfw;
m18_l_bwfw_z= l_bwfw_z;
m18_e_bwfw  = e_bwfw;
m18_e_bwfw_z= e_bwfw_z;
        
clearvars l_fwbw l_fwbw_z e_fwbw e_fwbw_z l_bwfw l_bwfw_z e_bwfw e_bwfw_z


