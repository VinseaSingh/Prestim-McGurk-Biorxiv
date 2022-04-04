% Code for Figure 3 is written by Vinsea AV Singh 
%
% Original paper:
%'Characterizing Variability of Audiovisual Speech Perception Based on 
% Periodic and Aperiodic Features of Pre-stimulus Brain Activity'
% Paper link - https://www.biorxiv.org/content/10.1101/2022.01.20.477172v1

%% This code plots figure 3a: Behavior (Inter-individual variability)

%{
   Input: 'McG_ta_PercentResponse.mat' - total percentage of McGurk /ta/ illusory response of all the participants
   Output: Scatter plot showing distribituion of total %age of McGurk illusion response for all participants
%}
clear; clc
load('McG_ta_PercentResponse.mat')

% Sort the percentage of /ta/ responses in ascending order
Sort_McG_ta = sort(ta);

% Find the participants with less than and greater than 50% Mcgurk perception 
below50 = find(Sort_McG_ta <= 50); % Rare group of perceivers
above50 = find(Sort_McG_ta >= 50); % Frequent group of perceivers

% Plot
x_rare = below50;                             %participant number 
y_rare = Sort_McG_ta(1:length(below50));      %percentage of McGurk percept /ta/

x_freq = above50;                             %participant number
y_freq = Sort_McG_ta(length(below50)+1:end);  %percentage of McGurk percept /ta/

figure;
scatter(x_rare,y_rare,80,'filled','o','r'); hold on;
scatter(x_freq,y_freq,80,'filled','o','b');
yline(50,'--', 'LineWidth',1.5)
title('(a) Inter-Individual Variability','FontSize',18);
xlabel('Participant number','FontSize',16); ylabel('% McGurk percept /ta/','FontSize',16);
legend('Rare perceivers','Frequent perceivers','FontSize',11); legend('boxoff')


%% This code plots figure 3b: Behavior (Inter-trial variability)

%{
   Input: 'McG_ta_PercentResponse.mat' - total percentage of McGurk /ta/ illusory response of all the participants
          'Subject_Percent_Cong.mat' - total percentage of congruent /pa/, /ta/, /ka/, and other 
                                       responses of all the participants during congruent trials
          'Subject_Percent_McG.mat' - total percentage of McGurk /pa/ (unisensory), /ta/ (illusory), 
                                      /ka/, and other responses of all the participants during McGurk trials

   Output: Violin plot showing inter-trial variability during incongruent and congruent stimulus 
           responses for both the rare and frequent group of perceivers

   Toolbox used for violin plot is 'Violin Plots for Matlab' by Bechtold and Bastian, 2016 
   (https://github.com/bastibe/Violinplot-Matlab) 
%}

% Load the required files in the folder "Data"
load('McG_ta_PercentResponse.mat')
load('Subject_Percent_Cong.mat') 
load('Subject_Percent_McG.mat')

% Add the violinplot toolbox to the MATLAB path
% All toolboxes are present in the folder "Toolboxes"
addpath('Violinplot-Matlab-master') 

% find index of rarely (<50%) and frequently (>50%) perceiving participants
rare_idx = find(ta <= 50);
freq_idx = find(ta > 50);

% Extract the percentage responses for the rare participants
for i = 1:length(rare_idx)
    sub = rare_idx(i);
    Cong_pa_rare(i) = Sub_Percent_Cong{sub}(1);
    Cong_ta_rare(i) = Sub_Percent_Cong{sub}(2);
    Cong_ka_rare(i) = Sub_Percent_Cong{sub}(3);
    McG_pa_rare(i)  = Sub_Percent_McG{sub}(1);
    McG_ta_rare(i)  = Sub_Percent_McG{sub}(2);
end

% Calculate mean response
m_Cong_pa_rare = mean(Cong_pa_rare);
m_Cong_ta_rare = mean(Cong_ta_rare);
m_Cong_ka_rare = mean(Cong_ka_rare);

m_rare.Cong = (m_Cong_pa_rare + m_Cong_ta_rare + m_Cong_ka_rare)./3;    
m_rare.McG_pa = mean(McG_pa_rare);
m_rare.McG_ta = mean(McG_ta_rare);

% Calculate standard deviation
std_pa_rare = std(Cong_pa_rare);
std_ta_rare = std(Cong_ta_rare);
std_ka_rare = std(Cong_ka_rare); 

sd_rare.Cong = (std_pa_rare + std_ta_rare + std_ka_rare)./3;                    
sd_rare.McG_pa = std(McG_pa_rare);            
sd_rare.McG_ta = std(McG_ta_rare);   


% Extract the percentage responses for the frequent participants
for i = 1:length(freq_idx)
    sub = freq_idx(i);
    Cong_pa_freq(i) = Sub_Percent_Cong{sub}(1);
    Cong_ta_freq(i) = Sub_Percent_Cong{sub}(2);
    Cong_ka_freq(i) = Sub_Percent_Cong{sub}(3);
    McG_pa_freq(i)  = Sub_Percent_McG{sub}(1);
    McG_ta_freq(i)  = Sub_Percent_McG{sub}(2);
end

% Calculate mean response
m_Cong_pa_freq = mean(Cong_pa_freq);
m_Cong_ta_freq = mean(Cong_ta_freq);
m_Cong_ka_freq = mean(Cong_ka_freq);

m_freq.Cong = (m_Cong_pa_freq + m_Cong_ta_freq + m_Cong_ka_freq)./3;    
m_freq.McG_pa = mean(McG_pa_freq);
m_freq.McG_ta = mean(McG_ta_freq);

% Calculate standard deviation
std_pa_freq = std(Cong_pa_freq);
std_ta_freq = std(Cong_ta_freq);
std_ka_freq = std(Cong_ka_freq); 

sd_freq.Cong = (std_pa_freq + std_ta_freq + std_ka_freq)./3;                    
sd_freq.McG_pa = std(McG_pa_freq);            
sd_freq.McG_ta = std(McG_ta_freq);  

% Plot
figure;
subplot(2,2,[1,2])
violinplot([McG_pa_rare',McG_ta_rare', Cong_pa_rare', Cong_ta_rare', Cong_ka_rare'],{'unisensory /pa/','Illusory /ta/','/pa/-/pa/','/ta/-/ta/','/ka/-/ka/'})
ylabel('Percentage responses','FontSize',14);
title('Rare perceivers','FontSize',16);

subplot(2,2,[3,4])
violinplot([McG_pa_freq',McG_ta_freq', Cong_pa_freq', Cong_ta_freq', Cong_ka_freq'],{'unisensory /pa/','Illusory /ta/','/pa/-/pa/','/ta/-/ta/','/ka/-/ka/'})
ylabel('Percentage responses','FontSize',14);
title('Frequent perceivers','FontSize',16);
sgtitle('(b) Inter-trial Variability','FontSize',18)


% significance testing of difference between percentage /ta/ and /pa/ percept during McGurk stimulus 
[h1, p1, ci1, stats1] = ttest2(McG_ta_rare, McG_pa_rare);
[h2, p2, ci2, stats2] = ttest2(McG_ta_freq, McG_pa_freq);

