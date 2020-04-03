    % author lfkrieger
% adjusted to experimental group CLm13-18
% start adjustments 04/2020

%% Part1
%% read-in .csv file
% from DeepLabCut output
% coordinates of mouse landmarks (less than with CLm1-12)
% 1-3 nose, 5-7 right ear, 8-10 left ear, 11-13 bodymass, 14-16 redLED_OFF,
% 17-19 blueLED_OFF, 20-22 redLED_ON, 23-25 blueLED_ON
file ='CLm15_20200220_CLm15_20200220_CLm15_20200220_CLm15_20200220_CLm15_20200220_Run6DeepCut_resnet50_CLmice13-18mar30shuffle1_350000.csv';
T = readtable(file);
T = table2array(T(3:end,:));

%% manually define TS identity accoring to order in time
% 0 = start/stop
% w/o laser
    % fwbw = 1; bwfw = 2
    % fwfw = -1; bwbw = -2 
% w laser
    % fwbw = 3; bwfw = 4
    % fwfw = -3; bwbw = -4 
% normally sin5: TSid = [0,1,2,1,2,1,2,1,0];
TSid = [0,1,2,1,2,1,4,1,2];

%% find direction TS
%red LED
    Tredon = str2double([T(:,22)]);
    Tredoff = str2double([T(:,16)]);

    for i = 1:length(Tredon)-1
        l(i) = (Tredon(i+1,1)-Tredon(i,1));    
    end
    RedOn = find(l>0.95); %frame on %+1 frame
    RedOff = find(l<-0.95); %frame off

%blue LED 
    Tblueon = str2double([T(:,25)]);
    Tblueoff = str2double([T(:,19)]);

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

% join red+blue LED TS
    dTS = [RedOn BlueOn];
    dTS = sort(dTS);
    dTS = vertcat(dTS,TSid);
    dTS = dTS';
    
%%
T = readtable(file); %because I want the original table be saved 
m15.D20200220.Run6={T,dTS};
clear ans BlueOn BlueOff dTS i RedOn RedOff T Tblueon Tblueoff Tredon Tredoff TSid file

% m15.D20200227.Run5{1}=T;  

%% Part2 
%% post-processing
%1.) imregister

img1 = imread('F:\CLm1-10_videos\CLm6\20191101\Run1\acA1920-25uc__21817890__20191101_142335857_0001.tiff');
img2 = imread('F:\CLm1-10_videos\CLm6\20191109\Run1\acA1920-25uc__21817890__20191109_134203058_0001.tiff');
img1 = rgb2gray(img1);
img2 = rgb2gray(img2);
imshowpair(img1, img2,'Scaling','joint')
[optimizer, metric] = imregconfig('multimodal') %or 'monomodal'
movingRegistered = imregister(img1, img2, 'affine', optimizer, metric);
[moving_reg,R_reg] = imregister(img1,img2,'affine',optimizer,metric);
figure
imshowpair(img1, movingRegistered,'Scaling','joint')

        %explanation of data structure:
        %strct containing table
        % to make calculation, convert to double: table2array (tabel to cell),
        % str2double (cell to double)
        % once done, convert back to original strucutre: num2cell(double to cell), then cell2table (cell to table) 

%2.) normalize coordinates to zero position = mouse sitting still before
%run onset or in first 150frames
fns = fieldnames(m6);
for i = 1:length(fns) % for all days
    fns2 = fieldnames(m6.(fns{i}));
    vec = [2,3,5,6,8,9,11,12,14,15,17,18,20,21,23,24,26,27,29,30,32,33];
    for j = 1:length(fns2) % for all Runs
        T = table2array(m6.(fns{i}).(fns2{j}){1,1}(1:end,:));
        for k = 1:length(vec)% for all x y columns
            %normalize to mean of first 150 frames
            %TS0 = m6.(fns{i}).(fns2{j}){1,2}(1,1); %alternatively all
            %frames before first TS
            M = mean(str2double(table2array(m6.(fns{i}).(fns2{j}){1,1}(3:152,vec(k)))));
            T(3:end,vec(k)) = num2cell(str2double(table2array(m6.(fns{i}).(fns2{j}){1,1}(3:end,vec(k))))-M);
        end
        m6.(fns{i}).(fns2{j}){1,1}=cell2table(T);
    end
end

%3.) excluding coord with low estimation-likelihood =NaN
% but what cutoff value? 95%
fns = fieldnames(m6);
for i = 1:length(fns) % for all days
    fns2 = fieldnames(m6.(fns{i}));
    for j = 1:length(fns2) % for all Runs
        T = (table2array(m6.(fns{i}).(fns2{j}){1,1}(1:end,:)));
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
        %dummy2=(str2double(table2array(m6.D20191107.Run4{1,1}(dummy,7))));
        %histogram(dummy2,50)
        m6.(fns{i}).(fns2{j}){1,1}=cell2table(T);
    end
end

%% get all TS of a type and associated timeseries of landmark coord.

%here for two landmarks, the center of bodymass and the right ear
l_fwbw =cell(4,2); %{1,1} = bodymass x, {1,2} = bodymass y
l_bwfw =cell(4,2); %{1,1} = bodymass x, {1,2} = bodymass y
e_fwbw =cell(4,2); %{1,1} = right ear x, {1,2} = right ear y
e_bwfw =cell(4,2); %{1,1} = right ear x, {1,2} = right ear y

fns = fieldnames(m6);
for i = 1:length(fns) % for all days
    fns2 = fieldnames(m6.(fns{i}));
    for j = 1:length(fns2) % for all Runs
    % find all TS of certain kind + associated frame number    
        %fwbw w/o laser
        idx_fw = find(m6.(fns{i}).(fns2{j}){1,2}(:,2) ==1);
        idx_fw2 = m6.(fns{i}).(fns2{j}){1,2}(idx_fw,1)+3;
        %bwfw w/o laser
        idx_bw = find(m6.(fns{i}).(fns2{j}){1,2}(:,2) ==2);
        idx_bw2 = m6.(fns{i}).(fns2{j}){1,2}(idx_bw,1)+3; 
        %fwbw w/ laser
        idx_fwL = find(m6.(fns{i}).(fns2{j}){1,2}(:,2) ==3);
        idx_fwL2 = m6.(fns{i}).(fns2{j}){1,2}(idx_fwL,1)+3;
        %bwfw w/laser
        idx_bwL = find(m6.(fns{i}).(fns2{j}){1,2}(:,2) ==4);
        idx_bwL2 = m6.(fns{i}).(fns2{j}){1,2}(idx_bwL,1)+3;
        %fwfw w/o laser
        idx_fwfw = find(m6.(fns{i}).(fns2{j}){1,2}(:,2) ==-1);
        idx_fwfw2 = m6.(fns{i}).(fns2{j}){1,2}(idx_fwfw,1)+3;
        %bwbw w/o laser
        idx_bwbw = find(m6.(fns{i}).(fns2{j}){1,2}(:,2) ==-2);
        idx_bwbw2 = m6.(fns{i}).(fns2{j}){1,2}(idx_bwbw,1)+3;
        %fwfw w/ laser
        idx_fwfwL = find(m6.(fns{i}).(fns2{j}){1,2}(:,2) ==-3);
        idx_fwfwL2 = m6.(fns{i}).(fns2{j}){1,2}(idx_fwfwL,1)+3;
        %bwbw w/ laser
        idx_bwbwL = find(m6.(fns{i}).(fns2{j}){1,2}(:,2) ==-4);
        idx_bwbwL2 = m6.(fns{i}).(fns2{j}){1,2}(idx_bwbwL,1)+3;
    % find tracking coordinates of body mass and right ear around TS
        T = table2array(m6.(fns{i}).(fns2{j}){1,1}(1:end,:));
        %l_fwbw AND e_fwbw
        for k = 1:length(idx_fw2) %normal fwbw x and y coord 
            l_fwbw{1,1}(:,size(l_fwbw{1,1},2)+1)= (T(idx_fw2(k)-40*6:idx_fw2(k)+40*6,17)); % bodymass x
            l_fwbw{1,2}(:,size(l_fwbw{1,2},2)+1)= (T(idx_fw2(k)-40*6:idx_fw2(k)+40*6,18)); % bodymass y
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
clearvars -except m6 l_bwfw l_fwbw e_bwfw e_fwbw

%% plot landmark coordinates aligned to direction TS

% here has to be done for each landmark individually
% raw values
    var = e_fwbw;
    var2 = 'fwbw ear x';
    var3 = 'fwbw ear y';
        figure %fwbw
        subplot(6,2,1)%fwbw RR protocol
        x=0:480; y=linspace(16,18,80); y=[y linspace(5,5,160)]; y=[y linspace(-10,-16,241)]; y2=NaN(1,240); y2=[y2 linspace(10,16,241)]; 
        plot(x,y,'color', [0.75 0.75 0.75]), hold on, plot(x,y2,'r'), vline([80 241],{'b','b'}), hline(0,'k'), xlim([0 481]), set(gca,'XTick',[],'YTick',[-15 15]); ylabel('rpm');
        subplot(6,2,2)
        plot(x,y,'color', [0.75 0.75 0.75]), hold on, plot(x,y2,'r'),vline([80 241],{'b','b'}), hline(0,'k'), xlim([0 481]), set(gca,'XTick',[],'YTick',[-15 15]); ylabel('rpm');
        subplot(6,2,[3 5 7 9 11]) %x coord
        plot(cell2mat(var{1,1}),':', 'color', [0.75 0.75 0.75]), hold on
        plot(cell2mat(var{2,1}), 'r','LineWidth',1), hold on
        plot(cell2mat(var{3,1}), 'r','LineWidth',1), hold on
        plot(cell2mat(var{4,1}), 'r','LineWidth',1)
        plot(nanmean(cell2mat(var{1,1})'),'b','LineWidth',2);
        vline([80 241],{'b','b'}), xlim([0 481]), ylim([-40 40]), title(var2), ylabel('lateral coord.'), xlabel('[seconds]');
        ax = gca;    
        ax.XTick = [0 40 80 120 160 200 240 280 320 360 400 440 480];        
        ax.XTickLabel = {'-6','-5','-4','-3','-2','-1','0','1','2','3','4','5','6'};       
        subplot(6,2,[4 6 8 10 12]) % y coord
        plot(cell2mat(var{1,2}),':', 'color', [0.75 0.75 0.75]), hold on
        plot(cell2mat(var{2,2}), 'r','LineWidth',1), hold on
        plot(cell2mat(var{3,2}), 'r','LineWidth',1), hold on
        plot(cell2mat(var{4,2}), 'r','LineWidth',1)
        plot(nanmean(cell2mat(var{1,2})'),'b','LineWidth',2);
        vline([80 241],{'b','b'}), xlim([0 481]), ylim([-150 100]) , title(var3), ylabel('longitudinal coord.'), xlabel('[seconds]');
        ax = gca;    
        ax.XTick = [0 40 80 120 160 200 240 280 320 360 400 440 480];        
        ax.XTickLabel = {'-6','-5','-4','-3','-2','-1','0','1','2','3','4','5','6'};
        
 
        
    var = e_bwfw;
    var2 = 'bwfw ear x';
    var3 = 'bwfw ear y';
        figure %bwfw
        subplot(6,2,1)%bwfw RR protocol
        x=0:480; y=linspace(-16,-18,80); y=[y linspace(-5,-5,160)]; y=[y linspace(10,16,241)]; y2=NaN(1,240); y2=[y2 linspace(-10,-16,241)]; 
        plot(x,y), hold on, plot(x,y2,'r'), vline([80 241],{'b','b'}), hline(0,'k'), xlim([0 481]), set(gca,'XTick',[],'YTick',[-15 15]); ylabel('rpm');
        subplot(6,2,2)
        plot(x,y), hold on, plot(x,y2,'r'), vline([80 241],{'b','b'}), hline(0,'k'), xlim([0 481]), set(gca,'XTick',[],'YTick',[-15 15]); ylabel('rpm');
        subplot(6,2,[3 5 7 9 11]) %x coord
        plot(cell2mat(var{1,1}),':', 'color', [0.75 0.75 0.75]), hold on
        plot(cell2mat(var{2,1}), 'r','LineWidth',1), hold on
        plot(cell2mat(var{3,1}), 'r','LineWidth',1), hold on
        plot(cell2mat(var{4,1}), 'r','LineWidth',1)
        plot(nanmean(cell2mat(var{1,1})'),'b','LineWidth',2);
        xlim([0 481]),ylim([-40 40]), vline([80 241],{'b','b'}), title(var2), ylabel('lateral coord.'), xlabel('[seconds]');
        ax = gca;    
        ax.XTick = [0 40 80 120 160 200 240 280 320 360 400 440 480];        
        ax.XTickLabel = {'-6','-5','-4','-3','-2','-1','0','1','2','3','4','5','6'};
        subplot(6,2,[4 6 8 10 12]) %y coord
        plot(cell2mat(var{1,2}),':', 'color', [0.75 0.75 0.75]), hold on
        plot(cell2mat(var{2,2}), 'r','LineWidth',1), hold on
        plot(cell2mat(var{3,2}), 'r','LineWidth',1), hold on
        plot(cell2mat(var{4,2}), 'r','LineWidth',1)
        plot(nanmean(cell2mat(var{1,2})'),'b','LineWidth',2);
        xlim([0 481]), ylim([-100 81]), vline([80 241],{'b','b'}), title(var3), ylabel('longitudinal coord.'), xlabel('[seconds]');
        ax = gca;    
        ax.XTick = [0 40 80 120 160 200 240 280 320 360 400 440 480];        
        ax.XTickLabel = {'-6','-5','-4','-3','-2','-1','0','1','2','3','4','5','6'};
        
        
%% z-transformation and plot

% ztransformation of the whole cell array _fwbw or _bwfw
 var = l_fwbw; %! RENAME var_z output variable!   
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

% plot scatter ztransformed coord of repetitions per timepoint 
% !must change VAR: landmark and fwbw!
% why not as line plot?
    var = e_fwbw_z;
    var2 = 'fwbw ear x';
    var3 = 'fwbw ear y';
        figure %fwbw
        subplot(6,2,1)%fwbe RR protocol
        x=0:480; y=linspace(16,18,80); y=[y linspace(5,5,160)]; y=[y linspace(-10,-16,241)]; y2=NaN(1,240); y2=[y2 linspace(10,16,241)]; 
        plot(x,y,'color', [0.75 0.75 0.75]), hold on, plot(x,y2,'r'), vline([80 241],{'b','b'}), hline(0,'k'), xlim([0 481]), set(gca,'XTick',[],'YTick',[-15 15]); ylabel('rpm');
        subplot(6,2,2)
        plot(x,y,'color', [0.75 0.75 0.75]), hold on, plot(x,y2,'r'),vline([80 241],{'b','b'}), hline(0,'k'), xlim([0 481]), set(gca,'XTick',[],'YTick',[-15 15]); ylabel('rpm');
        subplot(6,2,[3 5 7 9 11]) %x coord
        [rows cols] = size(var{1,1});
        plot((var{1,1}),':', 'color', [0.75 0.75 0.75]), hold on
        plot((var{2,1}), 'r','LineWidth',1), hold on
        plot((var{3,1}), 'r','LineWidth',1), hold on
        plot((var{4,1}), 'r','LineWidth',1)
        line([0 rows],[2 2],'Color',[0 0 0])
        line([0 rows],[-2 -2],'Color',[0 0 0])
        vline([80 241],{'b','b'}), xlim([0 481]), title(var2), ylabel('ztransformed lateral coord.'), xlabel('[seconds]');
        ax = gca;    
        ax.XTick = [0 40 80 120 160 200 240 280 320 360 400 440 480];        
        ax.XTickLabel = {'-6','-5','-4','-3','-2','-1','0','1','2','3','4','5','6'};
        subplot(6,2,[4 6 8 10 12]) % y coord
        [rows cols] = size(var{1,2});
        plot((var{1,2}),':', 'color', [0.75 0.75 0.75]), hold on
        plot((var{2,2}), 'r','LineWidth',1), hold on
        plot((var{3,2}), 'r','LineWidth',1), hold on
        plot((var{4,2}), 'r','LineWidth',1)
        line([0 rows],[2 2],'Color',[0 0 0])
        line([0 rows],[-2 -2],'Color',[0 0 0])
        xlim([0 481]), ylim([-5 5]), vline([80 241],{'b','b'}), title(var3), ylabel('ztransformed longitudinal coord.'), xlabel('[seconds]');
        ax = gca;    
        ax.XTick = [0 40 80 120 160 200 240 280 320 360 400 440 480];        
        ax.XTickLabel = {'-6','-5','-4','-3','-2','-1','0','1','2','3','4','5','6'};

    var = e_bwfw_z;
    var2 = 'bwfw ear x';
    var3 = 'bwfw ear y';
        figure %bwfw
        subplot(6,2,1)%bwfw RR protocol
        x=0:480; y=linspace(-16,-18,80); y=[y linspace(-5,-5,160)]; y=[y linspace(10,16,241)]; y2=NaN(1,240); y2=[y2 linspace(-10,-16,241)]; 
        plot(x,y), hold on, plot(x,y2,'r'), vline([80 241],{'b','b'}), hline(0,'k'), xlim([0 481]), set(gca,'XTick',[],'YTick',[-15 15]); ylabel('rpm');
        subplot(6,2,2)
        plot(x,y), hold on, plot(x,y2,'r'), vline([80 241],{'b','b'}), hline(0,'k'), xlim([0 481]), set(gca,'XTick',[],'YTick',[-15 15]); ylabel('rpm');
        subplot(6,2,[3 5 7 9 11]) %x coord
        [rows cols] = size(var{1,1});
        plot((var{1,1}),':', 'color', [0.75 0.75 0.75]), hold on
        plot((var{2,1}), 'r','LineWidth',1), hold on
        plot((var{3,1}), 'r','LineWidth',1), hold on
        plot((var{4,1}), 'r','LineWidth',1)
        line([0 rows],[2 2],'Color',[0 0 0])
        line([0 rows],[-2 -2],'Color',[0 0 0])
        xlim([0 481]), ylim([-5 5]), vline([80 241],{'b','b'}), title(var2), ylabel('ztransformed lateral coord.'), xlabel('[seconds]');
        ax = gca;    
        ax.XTick = [0 40 80 120 160 200 240 280 320 360 400 440 480];        
        ax.XTickLabel = {'-6','-5','-4','-3','-2','-1','0','1','2','3','4','5','6'};
        subplot(6,2,[4 6 8 10 12]) % y coord
        [rows cols] = size(var{1,2});
        plot((var{1,2}),':', 'color', [0.75 0.75 0.75]), hold on
        plot((var{2,2}), 'r','LineWidth',1), hold on
        plot((var{3,2}), 'r','LineWidth',1), hold on
        plot((var{4,2}), 'r','LineWidth',1)
        line([0 rows],[2 2],'Color',[0 0 0])
        line([0 rows],[-2 -2],'Color',[0 0 0])
        vline([80 241],{'b','b'}), xlim([0 481]), title(var3), ylabel('ztransformed longitudinal coord.'), xlabel('[seconds]'); 
        ax = gca;    
        ax.XTick = [0 40 80 120 160 200 240 280 320 360 400 440 480];        
        ax.XTickLabel = {'-6','-5','-4','-3','-2','-1','0','1','2','3','4','5','6'};
  
        
%%    
%OLD: as scatter plot        
% plot scatter ztransformed coord of repetitions per timepoint 
    %!must change landmark and fwbw!
    %fwbw
    var = l_fwbw_z;
    var2 = 'fwbw bodymass x';
    var3 = 'fwbw bodymass y';
    figure 
    % landmark x-coordinates
        subplot(1,2,1) 
        [rows cols] = size(var{1,1});
        for i=1:rows
        scatter(repmat(i,[cols 1]),var{1,1}(i,:), [],[0.75 0.75 0.75])
        hold on
        end
        line([0 rows],[2 2],'Color',[0 0 0])
        line([0 rows],[-2 -2],'Color',[0 0 0])
        yl = ylim; 
        line([rows/2 rows/2],[yl(1) yl(2)],'Color',[0 0 0],'LineStyle','--')
        % not normalTS 
        for j=2:4
            [rows cols] = size(var{j,1});
            for i=1:rows
            scatter(repmat(i,[cols 1]),var{j,1}(i,:), [],'r')
            hold on
            end
        end
        title(var2), ylabel('ztransf x-coord'), xlabel('40 frames/s');
    % landmark y-coordinates    
        subplot(1,2,2) 
        [rows cols] = size(var{1,2});
        for i=1:rows
        scatter(repmat(i,[cols 1]),var{1,2}(i,:), [],[0.75 0.75 0.75])
        hold on
        end
        line([0 rows],[2 2],'Color',[0 0 0])
        line([0 rows],[-2 -2],'Color',[0 0 0])
        yl = ylim; 
        line([rows/2 rows/2],[yl(1) yl(2)],'Color',[0 0 0],'LineStyle','--')
        % not normalTS 
        for j=2:4
            [rows cols] = size(var{j,2});
            for i=1:rows
            scatter(repmat(i,[cols 1]),var{j,2}(i,:), [],'r')
            hold on
            end
        end
        title(var3), ylabel('z-transf y-coord'), xlabel('40 frames/s');
        hold off

    %bwfw
    var = l_bwfw_z;   
    var2 = 'bwfw bodymass x';
    var3 = 'bwfw bodymass y';
    figure
    % landmark x-coordinates
        subplot(1,2,1) 
        [rows cols] = size(var{1,1});
        for i=1:rows
        scatter(repmat(i,[cols 1]),var{1,1}(i,:), [],[0.75 0.75 0.75])
        hold on
        end
        line([0 rows],[2 2],'Color',[0 0 0])
        line([0 rows],[-2 -2],'Color',[0 0 0])
        yl = ylim; 
        line([rows/2 rows/2],[yl(1) yl(2)],'Color',[0 0 0],'LineStyle','--')
        % not normalTS 
        for j=2:4
            [rows cols] = size(var{j,1});
            for i=1:rows
            scatter(repmat(i,[cols 1]),var{j,1}(i,:), [],'r')
            hold on
            end
        end
        title(var2), ylabel('z-transf x-coord'), xlabel('40 frames/s');
        subplot(1,2,2) 
        
        % landmark y-coordinates           
        [rows cols] = size(var{1,2});
        for i=1:rows
        scatter(repmat(i,[cols 1]),var{1,2}(i,:), [],[0.75 0.75 0.75])
        hold on
        end
        line([0 rows],[2 2],'Color',[0 0 0])
        line([0 rows],[-2 -2],'Color',[0 0 0])
        yl = ylim; 
        line([rows/2 rows/2],[yl(1) yl(2)],'Color',[0 0 0],'LineStyle','--')
        % not normalTS 
        for j=2:4
            [rows cols] = size(var{j,2});
            for i=1:rows
            scatter(repmat(i,[cols 1]),var{j,2}(i,:), [],'r')
            hold on
            end
        end
        title(var3), ylabel('z-transf y-coord'), xlabel('40 frames/s');
        hold off
        
        