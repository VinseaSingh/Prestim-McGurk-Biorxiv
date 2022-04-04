% Code for Figure 7(a) and Supplementary Figure 6 is written by Vinsea AV Singh 
%
% Original paper:
%'Characterizing Variability of Audiovisual Speech Perception Based on 
% Periodic and Aperiodic Features of Pre-stimulus Brain Activity'
% Paper link - https://www.biorxiv.org/content/10.1101/2022.01.20.477172v1

%% This code plots figure 7a and Supplementary Figure 6: Global Coherence (Inter-individual variability)

%{
   Input: 'Prestim_data.m'  - pre-processed Prestimulus EEG data.
          'Poststim_data.m' - pre-processes Poststimulus EEG data.
          
   Output: 'GlobCoh_InterSub_pre' - Global coherence capturing group-wise variability in the prestimulus duration.
           'GlobCoh_InterSub_post' - Global coherence capturing group-wise variability in the post-stimulus duration.
           Line plots showing the average of global coherence distribution for rare and frequent group of perceivers, 
           for incongruent and congruent conditions.

   Toolbox used for estimating global coherence - 'Chronux: a platform for analyzing neural signals (Bokil H, et al.,2010)'
                                                   (http://chronux.org/)
   Function file required -> 'Function_TmAvgGlobCoh.m' - for computing global coherence at inter-individual level of variability
                          
%}
clear; clc
restoredefaultpath
addpath(genpath('chronux_2_12'))

load('Prestim_PowerSpec_data.mat')
load('Poststim_PowerSpec_data.mat')

%%%%--------------------------Global Coherence (Inter-individual variability)--------------------------
% parameters
win = 0.8;   % 0.8(for 800 msec) 
params.Fs = 1000;
params.tapers=[3 5]; 
params.pad = 0;
params.fpass = [0.1 45]; 

cfg.choice = 'InterSub';
cfg.interval = 50;

% pre-stimulus 
GlobCoh_InterSub_pre = Function_TmAvgGlobCoh(Prestim_data, win, params, cfg);

% post-stimulus
GlobCoh_InterSub_post = Function_TmAvgGlobCoh(Poststim_data, win, params, cfg);

%%%%%%%%%%%%%Plot%%%%%%%%%%%%%
%%%%(Figure 7a) McGurk 
% Prestimulus
figure;
subplot(1,2,1)
plot(GlobCoh_InterSub_pre.f, GlobCoh_InterSub_pre.Group.McG.Ctot{1},'r','LineWidth',1.5); hold on;
plot(GlobCoh_InterSub_pre.f, GlobCoh_InterSub_pre.Group.McG.Ctot{2},'b','LineWidth',1.5);
h = gca; h.Box = 'off'; xlabel('Frequency (Hz)'); ylabel('Global Coherence')
title('Prestimulus: McGurk'); 

% Poststimulus
subplot(1,2,2)
plot(GlobCoh_InterSub_post.f, GlobCoh_InterSub_post.Group.McG.Ctot{1},'r','LineWidth',1.5); hold on;
plot(GlobCoh_InterSub_post.f, GlobCoh_InterSub_post.Group.McG.Ctot{2},'b','LineWidth',1.5);
ylim([0.1 0.65]);h = gca; h.Box = 'off'; xlabel('Frequency (Hz)'); ylabel('Global Coherence')
title('Poststimulus: McGurk'); legend('<50%- Rare Perceivers','>50%- Frequent Perceivers');legend('boxoff')
sgtitle('Figure 7(a): Global coherence (Inter-individual Variability)')


%%%%(Supplementary Figure 6)
%%%Prestimulus
% Cong /pa/ 
figure;
subplot(3,2,1)
plot(GlobCoh_InterSub_pre.f, GlobCoh_InterSub_pre.Group.Cong_pa.Ctot{1},'r','LineWidth',1.5); hold on;
plot(GlobCoh_InterSub_pre.f, GlobCoh_InterSub_pre.Group.Cong_pa.Ctot{2},'b','LineWidth',1.5);
ylim([0.1 0.65]); h = gca; h.Box = 'off';xlabel('Frequency (Hz)'); ylabel('Global Coherence')
title('Congruent /pa/ (Prestimulus)'); 

% Cong /ta/ 
subplot(3,2,3)
plot(GlobCoh_InterSub_pre.f, GlobCoh_InterSub_pre.Group.Cong_ta.Ctot{1},'r','LineWidth',1.5); hold on;
plot(GlobCoh_InterSub_pre.f, GlobCoh_InterSub_pre.Group.Cong_ta.Ctot{2},'b','LineWidth',1.5);
ylim([0.1 0.65]); h = gca; h.Box = 'off'; xlabel('Frequency (Hz)'); ylabel('Global Coherence')
title('Congruent /ta/ (Prestimulus)'); 

% Cong /ka/ 
subplot(3,2,5)
plot(GlobCoh_InterSub_pre.f, GlobCoh_InterSub_pre.Group.Cong_ka.Ctot{1},'r','LineWidth',1.5); hold on;
plot(GlobCoh_InterSub_pre.f, GlobCoh_InterSub_pre.Group.Cong_ka.Ctot{2},'b','LineWidth',1.5);
ylim([0.1 0.65]); h = gca; h.Box = 'off'; xlabel('Frequency (Hz)'); ylabel('Global Coherence')
title('Congruent /ka/ (Prestimulus)'); 


%%%Poststimulus
% Cong /pa/ 
subplot(3,2,2)
plot(GlobCoh_InterSub_post.f, GlobCoh_InterSub_post.Group.Cong_pa.Ctot{1},'r','LineWidth',1.5); hold on;
plot(GlobCoh_InterSub_post.f, GlobCoh_InterSub_post.Group.Cong_pa.Ctot{2},'b','LineWidth',1.5);
ylim([0.1 0.65]); h = gca; h.Box = 'off';xlabel('Frequency (Hz)'); ylabel('Global Coherence')
title('Congruent /pa/ (Poststimulus)'); 

% Cong /ta/ 
subplot(3,2,4)
plot(GlobCoh_InterSub_post.f, GlobCoh_InterSub_post.Group.Cong_ta.Ctot{1},'r','LineWidth',1.5); hold on;
plot(GlobCoh_InterSub_post.f, GlobCoh_InterSub_post.Group.Cong_ta.Ctot{2},'b','LineWidth',1.5);
ylim([0.1 0.65]); h = gca; h.Box = 'off'; xlabel('Frequency (Hz)'); ylabel('Global Coherence')
title('Congruent /ta/ (Poststimulus)'); 

% Cong /ka/ 
subplot(3,2,6)
plot(GlobCoh_InterSub_post.f, GlobCoh_InterSub_post.Group.Cong_ka.Ctot{1},'r','LineWidth',1.5); hold on;
plot(GlobCoh_InterSub_post.f, GlobCoh_InterSub_post.Group.Cong_ka.Ctot{2},'b','LineWidth',1.5);
ylim([0.1 0.65]); h = gca; h.Box = 'off'; xlabel('Frequency (Hz)'); ylabel('Global Coherence')
title('Congruent /ka/ (Poststimulus)'); 
sgtitle('Global Coherence: Inter-Individual Variability')




