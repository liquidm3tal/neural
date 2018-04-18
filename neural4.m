clc
close all
clear all
n_subj = 1;      %change this to 4 if you want to run each subject from its folder giving 20 brain graphs (select all 9 runs each time it asks)
%global allsubs; allsubs = 0;     %change this to 1 ONLY if running all 36 runs from SAME folder, needs n_subj=1 as well
%choose_subject=1; %doesnt do anything yet. YET
ch=[11 51 56 60];   %if you want the second figure to only plot one or more channel, just type in in here

 %   area = ['C9' ; 'P9' ; 'PO' ; 'PO']; %Type 9 after the letter if there is no second letter for that channel
%subscript =['Z' 'Z'  '7'  '8'];

%for subj=1:n_subj
            clear signal state params
 [signal, state, params] = load_data;
% first
clc
        clear trial_start trial_target trial_result
n_channels = size(signal,2); 



trial_start=find(diff([state.Feedback',NaN])'~=0)+1;
trial_target=state.TargetCode(trial_start);     trial_target=double(trial_target);

trial_end=find(diff([state.ResultCode',NaN])'~=0)+1; 
trial_result=state.ResultCode(trial_end);    trial_result=double(trial_result);

discard_ind=abs(trial_result-trial_target);
n_trials=length(trial_target); n_correct_trials=n_trials-sum(discard_ind);
accuracy_percentage = 100*(n_trials-sum(discard_ind))/n_trials;
correct_trial_start = trial_start(discard_ind~=1);  
correct_trial_end = trial_end(discard_ind~=1);
correct_trial_target = trial_target(discard_ind~=1);

% CT1_start = correct_trial_start(correct_trial_target==1);
% CT1_end = correct_trial_end(correct_trial_target==1);
% CT2_start = correct_trial_start(correct_trial_target==2);
% CT2_end = correct_trial_end(correct_trial_target==2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% wintool() may be useful here
win=400; %milliseconds
step=50; %ms
Fs=params.SamplingRate.NumericValue; %sampling freq in hertz
dt = 1/Fs; %sampling period
n_samp_per_trial=mean(correct_trial_end-correct_trial_start);
N=n_samp_per_trial; %length of sample signal data
t=[0:N-1]*dt; %time vector
t_total_trial = n_samp_per_trial/160; %in seconds
n_windows = floor(( n_samp_per_trial-win )/step);

for j=1:n_correct_trials
for i =1:n_windows
    begin = correct_trial_start(j)+(i-1)*step;
    finish= correct_trial_start(j)+(i-1)*step+win;
    ffts_to_avg(:,:,i) = abs(  fft( signal(begin:finish-1,:),Fs )  );
end
fft_trial_avg(:,:,j)=mean(ffts_to_avg,3)';
end

fft_single = fft_trial_avg(:,1:Fs/2,:);
fft_single_tar1 = fft_single(:,:,correct_trial_target==1);
fft_single_tar2 = fft_single(:,:,correct_trial_target==2);

[r2, amp1, amp2] = calc_rsqu(fft_single_tar1, fft_single_tar2);
figure(1)
surf([1:80],[1:64],r2);view(2);box on;axis tight;
colorbar; caxis([0 max(max(r2))]);
xlabel('f')
ylabel('channel')



%% Topoplot
figure(2)
r2_24 = r2(:,1:24);
r2_24 = mean(r2_24,2);
global amin amax; amin=0; amax=max(max(r2));
topoplot_mm(r2_24(:,:),'eloc64.txt','eeg')

































