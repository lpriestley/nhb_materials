%% Timecourse fMRI analysis %%
clearvars -except action_peaks offer_peaks; close all; clc;

base_path = ''; % Define base path to GitHub folder
cd(fullfile(base_path, 'nat_behav_submission/data/behavioural/'))
addpath(fullfile(base_path, 'fmt'))

GLM = '4.1';
stat = 1;
textout = 0;
getpeak = 0; 

timelock = 'offer-timelocked';
if strcmp(timelock, 'offer-timelocked')
    epochSuffix = '_epoched_offer';
elseif strcmp(timelock, 'action-timelocked')
    epochSuffix = '_epoched_action';
end

timecourse_dir = fullfile(base_path, 'nat_behav_submission/data/fMRI/timecourse/');
tc_graph_dir = fullfile(base_path, 'nat_behav_submission/data/fMRI/timecourse_outputs/');
feat_name = 'GLM01';

write_data = 0;


%%
subjects = {'S103','S104','S105','S106','S107','S109','S110','S111','S112', 'S114', 'S115','S116','S117','S118','S119','S120','S121','S122','S123','S124','S125','S126', 'S128', 'S130', 'S131','S132','S133','S135'};

if ismember(GLM, {'4.2a', '4.2b', '4.2c', '4.2d', '4.2a_low', '4.2c_low'})
    roi = {'DRN_custom'}; 
end

if ismember(GLM, {'4.1', '4.1b','S6'})
    roi = {'DRN_custom', 'MBD', 'HB', 'BF', 'LC_Pauli'}; 
end

if ismember(GLM, {'S3a', 'S3b'})
    roi = {'MBD'}; 
end

if ismember(GLM, {'4.3'})
    seed_roi = {'dACCsphere7'}; roi = {'DRN_custom'}; 
end

behavAll = csvread('behavioural_data_for_tc.csv', 1, 0);

%%
for r = 1:numel(roi)
    
    z = 0;
    trial_data = [];
    
    %% get the action and stimulus onset time from behavioural data
    
    subjects = unique(behavAll(:,1)); % List of subjects
    
    for is=1:length(subjects)
        
        z = z + 1;
        s = subjects(is);
        
        % load data
        subStr = ['S', num2str(s)];
        
        % prepare behaviour
        dataBehav = behavAll(behavAll(:,1) == s, :);
        
        % remove trials in subjects that scanner crashed before the end
        if      s == 103 %remove the last trial
            dataBehav(end,:) = [];
        elseif  s == 104 %remove the last 7 trials
            dataBehav((end-6:end),:) = [];
        end
        
        trial = dataBehav(:,2);
        offer = dataBehav(:,3);
        decision = dataBehav(:,4);
        mu_val = dataBehav(:, 13);
        env_bin = dataBehav(:,28);
        prev_action = dataBehav(:, 46);
        prev_policy = dataBehav(:, 47);
        mu_val_pe = offer - mu_val;
        
        policy_change = NaN(1, length(trial))';
        for i=1:length(policy_change)
            if ~isnan(prev_policy(i))
                if prev_policy(i)==decision(i)
                    policy_change(i)=0;
                elseif prev_policy(i)~=decision(i)
                    policy_change(i)=1;
                end
            end
        end

        action_change = NaN(1, length(trial))';
        for i=1:length(action_change)
            if ~isnan(prev_action(i))
                if prev_action(i)==decision(i)
                    action_change(i)=0;
                elseif prev_action(i)~=decision(i)
                    action_change(i)=1;
                end
            end
        end

        pursue_change = NaN(1, length(trial))';
        for i=1:length(pursue_change)
            if ~isnan(prev_policy(i))
                if prev_policy(i)==0 & decision(i)==1
                    pursue_change(i)=1;
                else
                    pursue_change(i)=0;
                end
            end
        end
        
        reject_change = NaN(1, length(trial))';
        for i=1:length(reject_change)
            if ~isnan(prev_policy(i))
                if prev_policy(i)==1 & decision(i)==0
                    reject_change(i)=1;
                else
                    reject_change(i)=0;
                end
            end
        end

        congruent_change = NaN(1, length(trial))';
        for i=1:length(congruent_change)
            if ~isnan(prev_policy(i))
                if reject_change(i)==1 & env_bin(i)==1 % if they switch to reject in the rich environment
                    congruent_change(i)=1;
                elseif pursue_change(i)==1 & env_bin(i)==-1% if they switch to pursue in the poor environment
                    congruent_change(i)=1;
                else
                    congruent_change(i)=0;
                end
            end
        end
        
        incongruent_change = NaN(1, length(trial))';
        for i=1:length(incongruent_change)
            if ~isnan(prev_policy(i))
                if reject_change(i)==1 & env_bin(i)==-1 % if they switch to reject in the poor environment
                    incongruent_change(i)=1;
                elseif pursue_change(i)==1 & env_bin(i)==1% if they switch to pursue in the rich environment
                    incongruent_change(i)=1;
                else
                    incongruent_change(i)=0;
                end
            end
        end
        
        congruent_vs_incongruent = NaN(1, length(trial))';
        for i=1:length(congruent_vs_incongruent)
            if ~isnan(prev_policy(i))
                if congruent_change(i)==1
                    congruent_vs_incongruent(i)=1;
                elseif incongruent_change(i)==1
                    congruent_vs_incongruent(i)=-1;
                else
                    congruent_vs_incongruent(i)=0;
                end
            end
        end
        
        congruent_poor = NaN(1, length(trial))';
        for i=1:length(congruent_poor)
            if ~isnan(prev_policy(i))
                if congruent_change(i)==1 & offer(i)==10 & env_bin(i)==-1
                    congruent_poor(i)=1;
                else
                    congruent_poor(i)=0;
                end
            end
        end
        
        congruent_poor_low = NaN(1, length(trial))';
        for i=1:length(congruent_poor_low)
            if ~isnan(prev_policy(i))
                if congruent_change(i)==1 & offer(i)==5 & env_bin(i)==-1
                    congruent_poor_low(i)=1;
                else
                    congruent_poor_low(i)=0;
                end
            end
        end

        incongruent_poor = NaN(1, length(trial))';
        for i=1:length(incongruent_poor)
            if ~isnan(prev_policy(i))
                if incongruent_change(i)==1 & offer(i)==10 & env_bin(i)==-1
                    incongruent_poor(i)=1;
                else
                    incongruent_poor(i)=0;
                end
            end
        end
        
        congruent_rich = NaN(1, length(trial))';
        for i=1:length(congruent_rich)
            if ~isnan(prev_policy(i))
                if congruent_change(i)==1 & offer(i)==10 & env_bin(i)==1
                    congruent_rich(i)=1;
                else
                    congruent_rich(i)=0;
                end
            end
        end
        
        congruent_rich_low = NaN(1, length(trial))';
        for i=1:length(congruent_rich_low)
            if ~isnan(prev_policy(i))
                if congruent_change(i)==1 & offer(i)==5 & env_bin(i)==1
                    congruent_rich_low(i)=1;
                else
                    congruent_rich_low(i)=0;
                end
            end
        end

        incongruent_rich = NaN(1, length(trial))';
        for i=1:length(incongruent_rich)
            if ~isnan(prev_policy(i))
                if incongruent_change(i)==1 & offer(i)==10 & env_bin(i)==1
                    incongruent_rich(i)=1;
                else
                    incongruent_rich(i)=0;
                end
            end
        end

        % Prepare regressors
        % main
        REG.trial = trial;
        REG.offer = offer;
        REG.decision = decision;
        REG.mu_val = mu_val;
        REG.env_bin = env_bin;
        REG.prev_action = prev_action;
        REG.prev_policy = prev_policy;
        
        REG.policy_change = policy_change;
        REG.action_change = action_change;
        
        REG.congruent_poor = congruent_poor; 
        REG.congruent_rich = congruent_rich; 
        REG.incongruent_poor = incongruent_poor; 
        REG.incongruent_rich = incongruent_rich; 
        REG.congruent_poor_low = congruent_poor_low; 
        REG.congruent_rich_low = congruent_rich_low; 

        REG.incongruent_change = incongruent_change; 
        REG.mu_val_pe = mu_val_pe;

        REG.constant = ones(length(REG.trial), 1);
        
        % load time-series data
        load([timecourse_dir,subStr,'/',feat_name,'/',(['tc_' roi{r} epochSuffix])]);
        
        
        %% Run GLMs
        
        if  strcmp(GLM, '4.1')
            
            % regressor design matrix (don't forget constant!)
            dmat =  [REG.policy_change, REG.trial, REG.constant];
            
            dmat(isnan(REG.policy_change),:)=[];
            trial_data(isnan(REG.policy_change),:)=[];
            
            % contrast design matrix
            contrasts = diag(ones(size(dmat,2),1));
            
            % normalise data
            dmat(:,1:(size(dmat,2)-1)) = zscore(dmat(:,1:(size(dmat,2)-1)));
            
            % beta X time output
            betaOut = ols(trial_data,dmat,contrasts);
            allBeta.(sprintf('%s',feat_name)).(sprintf('%s',roi{r}))(:,:,z) = betaOut;
            
            % write effect size
            if  textout
                filename = [timecourse_dir,s,'/',session,'/',feat_name,'/',([roi{r},'_effectSize'])];
                save(filename,'betaOut');
            end
        end

        if  strcmp(GLM, '4.1b')
            
            % regressor design matrix (don't forget constant!)
            dmat =  [REG.action_change, REG.trial, REG.constant];
            
            dmat(isnan(REG.action_change),:)=[];
            trial_data(isnan(REG.action_change),:)=[];
            
            % contrast design matrix
            contrasts = diag(ones(size(dmat,2),1));
            
            % normalise data
            dmat(:,1:(size(dmat,2)-1)) = zscore(dmat(:,1:(size(dmat,2)-1)));
            
            % beta X time output
            betaOut = ols(trial_data,dmat,contrasts);
            allBeta.(sprintf('%s',feat_name)).(sprintf('%s',roi{r}))(:,:,z) = betaOut;
            
            % write effect size
            if  textout
                filename = [timecourse_dir,s,'/',session,'/',feat_name,'/',([roi{r},'_effectSize'])];
                save(filename,'betaOut');
            end
        end
        
        if strcmp(GLM, '4.2a') % effect of congruent pursue change
            
            % regressor design matrix (don't forget constant!)
            dmat =  [REG.congruent_poor, REG.trial, REG.constant];
            
            dmat(isnan(REG.congruent_poor),:)=[];
            trial_data(isnan(REG.congruent_poor),:)=[];
            
            % filter mid-opt trials only
            mag_index = offer;
            mag_index(isnan(REG.congruent_poor),:)=[];
            
            % filter according to environment
            dmat = dmat(mag_index==10,:);
            trial_data = trial_data(mag_index==10,:);
            
            % contrast design matrix
            contrasts = diag(ones(size(dmat,2),1));
            
            % normalise data
            dmat(:,1:(size(dmat,2)-1)) = zscore(dmat(:,1:(size(dmat,2)-1)));
            
            % beta X time output
            betaOut = ols(trial_data,dmat,contrasts);
            allBeta.(sprintf('%s',feat_name)).(sprintf('%s',roi{r}))(:,:,z) = betaOut;
            
            % write effect size
            if  textout
                filename = [timecourse_dir,s,'/',session,'/',feat_name,'/',([roi{r},'_effectSize'])];
                save(filename,'betaOut');
            end
        end

        if  strcmp(GLM, '4.2b') % effect of incongruent pursue change

            % regressor design matrix (don't forget constant!)
            dmat =  [REG.incongruent_rich, REG.trial, REG.constant];

            dmat(isnan(REG.incongruent_rich),:)=[];
            trial_data(isnan(REG.incongruent_rich),:)=[];

            % contrast design matrix
            contrasts = diag(ones(size(dmat,2),1));

            % normalise data
            dmat(:,1:(size(dmat,2)-1)) = zscore(dmat(:,1:(size(dmat,2)-1)));

            % beta X time output
            betaOut = ols(trial_data,dmat,contrasts);
            allBeta.(sprintf('%s',feat_name)).(sprintf('%s',roi{r}))(:,:,z) = betaOut;

            % write effect size
            if  textout
                filename = [timecourse_dir,s,'/',session,'/',feat_name,'/',([roi{r},'_effectSize'])];
                save(filename,'betaOut');
            end
        end
        
        if  strcmp(GLM, '4.2c') % effect of congruent reject change
            
            % regressor design matrix (don't forget constant!)
            dmat =  [REG.congruent_rich, REG.trial, REG.constant];
            
            dmat(isnan(REG.congruent_rich),:)=[];
            trial_data(isnan(REG.congruent_rich),:)=[];
            
            % filter mid-opt trials only
            mag_index = offer;
            mag_index(isnan(REG.congruent_rich),:)=[];
            
            % filter according to environment
            dmat = dmat(mag_index==10,:);
            trial_data = trial_data(mag_index==10,:);
            
            % contrast design matrix
            contrasts = diag(ones(size(dmat,2),1));
            
            % normalise data
            dmat(:,1:(size(dmat,2)-1)) = zscore(dmat(:,1:(size(dmat,2)-1)));
            
            % beta X time output
            betaOut = ols(trial_data,dmat,contrasts);
            allBeta.(sprintf('%s',feat_name)).(sprintf('%s',roi{r}))(:,:,z) = betaOut;
            
            % write effect size
            if  textout
                filename = [timecourse_dir,s,'/',session,'/',feat_name,'/',([roi{r},'_effectSize'])];
                save(filename,'betaOut');
            end
        end

        if strcmp(GLM, '4.2d') % effect of incongruent reject change
            
            % regressor design matrix (don't forget constant!)
            dmat =  [REG.incongruent_poor, REG.trial, REG.constant];
            
            dmat(isnan(REG.incongruent_poor),:)=[];
            trial_data(isnan(REG.incongruent_poor),:)=[];
            
            % contrast design matrix
            contrasts = diag(ones(size(dmat,2),1));
            
            % normalise data
            dmat(:,1:(size(dmat,2)-1)) = zscore(dmat(:,1:(size(dmat,2)-1)));
            
            % beta X time output
            betaOut = ols(trial_data,dmat,contrasts);
            allBeta.(sprintf('%s',feat_name)).(sprintf('%s',roi{r}))(:,:,z) = betaOut;
            
            % write effect size
            if  textout
                filename = [timecourse_dir,s,'/',session,'/',feat_name,'/',([roi{r},'_effectSize'])];
                save(filename,'betaOut');
            end
        end
        
        if  strcmp(GLM, '4.2a_low') % effect of congruent reject change
            
            % regressor design matrix (don't forget constant!)
            dmat =  [REG.congruent_poor_low, REG.trial, REG.constant];
            
            dmat(isnan(REG.congruent_poor_low),:)=[];
            trial_data(isnan(REG.congruent_poor_low),:)=[];
            
            % contrast design matrix
            contrasts = diag(ones(size(dmat,2),1));
            
            % normalise data
            dmat(:,1:(size(dmat,2)-1)) = zscore(dmat(:,1:(size(dmat,2)-1)));
            
            % beta X time output
            betaOut = ols(trial_data,dmat,contrasts);
            allBeta.(sprintf('%s',feat_name)).(sprintf('%s',roi{r}))(:,:,z) = betaOut;
            
            % write effect size
            if  textout
                filename = [timecourse_dir,s,'/',session,'/',feat_name,'/',([roi{r},'_effectSize'])];
                save(filename,'betaOut');
            end
        end
        
        if strcmp(GLM, '4.2c_low') % effect of congruent reject change
            
            % regressor design matrix (don't forget constant!)
            dmat =  [REG.congruent_rich_low, REG.trial, REG.constant];
            
            dmat(isnan(REG.congruent_rich_low),:)=[];
            trial_data(isnan(REG.congruent_rich_low),:)=[];
            
            % contrast design matrix
            contrasts = diag(ones(size(dmat,2),1));
            
            % normalise data
            dmat(:,1:(size(dmat,2)-1)) = zscore(dmat(:,1:(size(dmat,2)-1)));
            
            % beta X time output
            betaOut = ols(trial_data,dmat,contrasts);
            allBeta.(sprintf('%s',feat_name)).(sprintf('%s',roi{r}))(:,:,z) = betaOut;
            
            % write effect size
            if  textout
                filename = [timecourse_dir,s,'/',session,'/',feat_name,'/',([roi{r},'_effectSize'])];
                save(filename,'betaOut');
            end
        end

        if  strcmp(GLM, '4.3') % PPI as a function of policy-switch
            
            % load time-series data
            seed_TC = load([timecourse_dir,subStr,'/',feat_name,'/',(['tc_' seed_roi{1} epochSuffix])]);
            roi_TC  = load([timecourse_dir,subStr,'/',feat_name,'/',(['tc_' roi{r} epochSuffix])]);
            seed_TC.trial_data(isnan(REG.env_bin),:)=[];
            
            % run GLM
            for i = 1:size(seed_TC.trial_data,2)
                
                % load ROI time-series data
                REG.TC  = roi_TC.trial_data(:,i);
                dmat = [REG.env_bin, REG.trial, REG.constant];
                dmat = [REG.TC, dmat];
                
                % remove trials with no response
                dmat(isnan(REG.env_bin),:)=[];
                
                % normalise data
                dmat(:,1:(size(dmat,2)-1)) = zscore(dmat(:,1:(size(dmat,2)-1)));
                
                % create PPI regressor
                REG.PPI = zscore (dmat(:,1) .* dmat(:,2));
                dmat = [REG.PPI, dmat];
                
                % contrast design matrix
                contrasts = diag(ones(size(dmat,2),1));
                
                % beta X time output
                betaOut(:,i) = ols(seed_TC.trial_data(:,i),dmat,contrasts);
                clear dmat contrasts REG.TC REG.PPI
            end
            allBeta.(sprintf('%s',feat_name)).(sprintf('%s',roi{r}))(:,:,z) = betaOut;
        end

        if  strcmp(GLM, 'S3a') % effect of response
            
            % regressor design matrix (don't forget constant!)
            dmat =  [REG.decision, REG.trial, REG.constant];

            % contrast design matrix
            contrasts = diag(ones(size(dmat,2),1));
            
            % normalise data
            dmat(:,1:(size(dmat,2)-1)) = zscore(dmat(:,1:(size(dmat,2)-1)));
            
            % beta X time output
            betaOut = ols(trial_data,dmat,contrasts);
            allBeta.(sprintf('%s',feat_name)).(sprintf('%s',roi{r}))(:,:,z) = betaOut;
            
            % write effect size
            if  textout
                filename = [timecourse_dir,s,'/',session,'/',feat_name,'/',([roi{r},'_effectSize'])];
                save(filename,'betaOut');
            end
        end
        
        if  strcmp(GLM, 'S3b') % effect of offer
            
            % regressor design matrix (don't forget constant!)
            dmat =  [REG.mu_val_pe, REG.decision, REG.trial, REG.constant];
            
            dmat(isnan(REG.mu_val_pe),:)=[];
            trial_data(isnan(REG.mu_val_pe),:)=[];
            
            % contrast design matrix
            contrasts = diag(ones(size(dmat,2),1));
            
            % normalise data
            dmat(:,1:(size(dmat,2)-1)) = zscore(dmat(:,1:(size(dmat,2)-1)));
            
            % beta X time output
            betaOut = ols(trial_data,dmat,contrasts);
            allBeta.(sprintf('%s',feat_name)).(sprintf('%s',roi{r}))(:,:,z) = betaOut;
            
            % write effect size
            if  textout
                filename = [timecourse_dir,s,'/',session,'/',feat_name,'/',([roi{r},'_effectSize'])];
                save(filename,'betaOut');
            end
        end

        if  strcmp (GLM, 'S6') 
            
            % regressor design matrix (don't forget constant!)
            dmat =  [REG.incongruent_change, REG.trial, REG.constant];
            
            dmat(isnan(REG.incongruent_change),:)=[];
            trial_data(isnan(REG.incongruent_change),:)=[];
            
            % contrast design matrix
            contrasts = diag(ones(size(dmat,2),1));
            
            % normalise data
            dmat(:,1:(size(dmat,2)-1)) = zscore(dmat(:,1:(size(dmat,2)-1)));
            
            % beta X time output
            betaOut = ols(trial_data,dmat,contrasts);
            allBeta.(sprintf('%s',feat_name)).(sprintf('%s',roi{r}))(:,:,z) = betaOut;
            
            % write effect size
            if  textout
                filename = [timecourse_dir,s,'/',session,'/',feat_name,'/',([roi{r},'_effectSize'])];
                save(filename,'betaOut');
            end
        end      
    end
end
%% write data out for graphs
contrast_to_plot = 1;

for r = 1:numel(roi)
    
    beta = squeeze(allBeta.(sprintf('%s',feat_name)).(sprintf('%s',roi{r}))(contrast_to_plot,:,:));
    
    if write_data
        tc_data = {};
        tc_data(:, 1) = num2cell(mean(beta'));
        tc_data(:, 2) = num2cell(std(beta')/sqrt(size(beta',1)));
        tc_data(:, 3) = roi(r);
        tc_data(:, 4) = {['GLM_0', num2str(GLM)]};
        tc_data(:, 5) = {['contr_', num2str(contrast_to_plot)]};
        if strcmp(GLM, '4.3') 
            filename = [tc_graph_dir, 'GLM_0', num2str(GLM),'_contr_', num2str(contrast_to_plot),'_',seed_roi{1},'_to_',roi{r},'_PPI_timecourse.csv'];
        else
            filename = [tc_graph_dir, 'GLM_0', num2str(GLM),'_contr_', num2str(contrast_to_plot),'_',roi{r},'_',timelock,'_timecourse.csv'];
        end
        writetable(cell2table(tc_data), filename);
        clear tc_data
    end
    
end


%% stats

if stat
    
    clearvars window LOOT peak
    numsession = 1:numel(subjects);
    pre_win = 2;
    post_win = 8;
    upsample = 20;
    TR = 1.775;
    nsamples = round(((pre_win+post_win)./TR)*upsample);
    start = 40;
    finish = max(nsamples);
    reg = contrast_to_plot;
    
    % find peak in a specified window
    for r = 1:numel(roi)
        
        for s = 1:length(numsession)
            
            numsession(s) = [];
            window  = mean(squeeze(allBeta.(sprintf('%s',feat_name)).(sprintf('%s',roi{r}))(reg,:,numsession)),2);
            
            [m,i]  = max(abs(window(start:finish)));
            LocPeak(s,r)=(i*10)/length(window);
            i      = i + (start-1);
            peak(s,r) = allBeta.(sprintf('%s',feat_name)).(sprintf('%s',roi{r}))(reg,i,s);
            numsession = 1:numel(subjects);
            clear window fw
            
        end
        [h,p(r),ci,stats]= ttest(peak(:,r));
        t(r) = stats.tstat;
    end

    [bonf_p, bonf_h] = bonf_holm(p)
    
    if write_data
        tb = num2cell(peak);
        tb(:, size(tb,2)+1) = {['GLM_0', num2str(GLM)]};
        tb(:, size(tb,2)+1) = {['contr', num2str(contrast_to_plot)]};
        tb(:, size(tb,2)+1) = {timelock};
        tb = cell2table(tb, 'VariableNames', [roi, 'GLM', 'contrast', 'timelock']);
        writetable(tb,[tc_graph_dir, 'GLM_0', num2str(GLM),'_contr_', num2str(contrast_to_plot),'_',timelock,'_peaks.csv']);
        clear tc_data
    end
end
