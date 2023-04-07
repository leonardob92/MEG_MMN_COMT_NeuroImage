%% COMT - PREDICTION ERROR - MEG
% NeuroImage paper:
%https://www.sciencedirect.com/science/article/pii/S1053811921002317

%%

%BEFORE PROCEEDING, PLEASE NOTE:

%#1
%Data can be provided only upon reasonable requests and when we are sure that is was completely anonymized and
%could not be linked to the original personal information of the participants.
%If interested, please contact me, Leonardo Bonetti: leonardo.bonetti(at).clin.au.dk

%#2
%This script reports the statistical analysis for MEG sensors and
%pre-processing and statistical analysis for MEG sources.
%The pre-processing at MEG sensor level followed a standard pipeline
%described in detail in the published paper.

%#3
%if you notice any typo or mistake, you are very welcome to let me know so
%that I can fix them.
%Of course, I am also available to provide more information about the codes.
%Thanks.

%#4
%Information on genotyping is not reported here, but descripted in details
%in the paper. The following codes relate to the MEG part only.

%% OSL - adding paths to functions

addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/osl/osl-core'); %this add the path (not necessary being in the osl-core directory if you have this one)
osl_startup
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/scripts_osl_learningbach'); %path to this specific folder (there are MEG_leonardo_functions for MEG_sensors)
% addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/Islands_new'); %path to islands and violin functions
% addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/Islands_new/schemaball-master'); %path to schemaball functions
% addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/altmany-export_fig-412662f'); %trying to get higher quality images

%% calculating t-tests (COMT - MMN)

%%% Data has been preprocessed using standard procedures, as described in the following paper:
%https://www.sciencedirect.com/science/article/pii/S1053811921002317
%Here, we report only the subsequent steps, which are the most meaningful
%to understand our results

subjs_list = dir('/scratch1/MINDLAB2015_MEG-TunteetMM/muMUFE_ns/Nader_avg_cmb/avg*');
subjs_homo = {'002','004','007','014','021','022','025','026','027','028','030','033','034','035','045','048','051','052','053','057','067','069','070','071','073','074','075','076','078','081','082','085','086','097','103','105','106','111','113','114','116','117','118','120','122','125','131','132','138','140'};
subjs_hetero = {'003','005','011','012','015','017','019','029','031','032','036','037','038','039','040','043','046','047','049','050','054','055','056','058','059','061','062','063','064','065','066','068','072','079','080','083','084','090','092','093','095','102','107','108','109','110','112','121','123','124','128','130','133','134','135','136','137','139'};

%reshaping data matching the original ID with the new sequential numbers
caz = string('E12');
subjs_homo2 = cell(length(subjs_homo),1);
for ii = 1:length(subjs_homo)
    subjs_homo2(ii) = {strcat(caz,string(subjs_homo{ii}))};
end
clear homo2
count = 0;
for ii = 1:length(subjs_homo2)
    for pp = 1:length(subjs_list)        
        if contains(string(subjs_list(pp).name),string(subjs_homo2{ii}))
            count = count + 1;
            homo2(count,1) = pp;
        end
    end
end
subjs_hetero2 = cell(length(subjs_hetero),1);
for ii = 1:length(subjs_hetero2)
    subjs_hetero2(ii) = {strcat(caz,string(subjs_hetero{ii}))};
end
clear hetero2
count = 0;
for ii = 1:length(subjs_hetero2)
    for pp = 1:length(subjs_list)        
        if contains(string(subjs_list(pp).name),string(subjs_hetero2{ii}))
            count = count + 1;
            hetero2(count,1) = pp;
        end
    end
end
%work around the data that has been saved with different numbers.. (TO BE RUN ONLY THE FIRST TIME!)
first_time = 0;
if first_time == 1
    for ii = 1:length(subjs_list)
        klm = subjs_list(ii).name;
        load(subjs_list(ii).name);
        save(['2' subjs_list(ii).name],'-struct',klm(1:end-4));
    end
end
subjs_list2 = dir('/scratch1/MINDLAB2015_MEG-TunteetMM/muMUFE_ns/Nader_avg_cmb/2avg*');

%% LOAD PROPER DATA AND THEN CALCULATE T-TESTS (INDEPENDENTLY FOR EACH DEVIANT OR AVERAGING THE DEVIANTS TOGETHER)
%LAST UPDATED 04/06/2019 

deviant = 5; %insert the number of deviant that you want to compute (insert 0 if you want to average all of the deviants together) (1 = intensity; 2 = location; 3 = pitch; 4 = rhythm; 5 = slide; 6 = timbre)
% timel = [31:122]; %time for MCS analysis
timel = [1:151]; %time for waveform plotting purposes

if deviant ~= 0 
    datahomo = zeros(204,length(timel),length(subjs_homo2));
    countd = deviant;
    for ii = 1:length(homo2)/6
        b = load([subjs_list2(homo2(countd)).folder '/' subjs_list2(homo2(countd)).name],'avg');
        datahomo(1:204,:,ii) = b.avg(1:204,timel);
        countd = countd + 6;
        disp(['just done subj number ' num2str(ii)])
    end
    datahetero = zeros(204,length(timel),length(subjs_hetero2));
    countd = deviant;
    for ii = 1:length(hetero2)/6
        b = load([subjs_list2(hetero2(countd)).folder '/' subjs_list2(hetero2(countd)).name],'avg');
        datahetero(1:204,:,ii) = b.avg(1:204,timel);
        countd = countd + 6;
        disp(['just done subj number ' num2str(ii)])
    end
else
    %hetero
    count = 0;
    dataeteroav = zeros(204,length(timel),length(subjs_hetero));
    for ii = 1:length(subjs_hetero)
        b = zeros(204,length(timel),6); %avg dimensions for each of the 6 deviants
        for jj = 1:6
            count = count + 1;
            a = load([subjs_list2(hetero2(count)).folder '/' subjs_list2(hetero2(count)).name],'avg');
            b(:,:,jj) = a.avg(1:204,timel);
        end
        dataeteroav(:,:,ii) = mean(b,3);
        disp(['just done subj number ' num2str(ii)])
    end
    %homo
    count = 0;
    datahomoav = zeros(204,length(timel),length(subjs_homo));
    for ii = 1:length(subjs_homo)
        b = zeros(204,length(timel),6); %avg dimensions for each of the 6 deviants
        for jj = 1:6
            count = count + 1;
            a = load([subjs_list2(homo2(count)).folder '/' subjs_list2(homo2(count)).name],'avg');
            b(:,:,jj) = a.avg(1:204,timel);
        end
        datahomoav(:,:,ii) = mean(b,3);
        disp(['just done subj number ' num2str(ii)])
    end
end
% actual two-sample t-tests
if deviant == 0
    SS = size(dataeteroav);
else
    SS = size(datahetero);
end
DATAP = zeros(SS(1),SS(2));
TSTAT = zeros(SS(1),SS(2));
for ii = 1:SS(1)
    for jj = 1:SS(2)
        if deviant == 0
            a = squeeze(dataeteroav(ii,jj,:));
            b = squeeze(datahomoav(ii,jj,:));
        else
            a = squeeze(datahetero(ii,jj,:));
            b = squeeze(datahomo(ii,jj,:));
        end
        [~,p,~,stats] = ttest2(a,b);
        DATAP(ii,jj) = p;
        TSTAT(ii,jj) = stats.tstat;
    end
end
% save averageddeviants_stats.mat DATAP TSTAT

%% working on data to make it suitable for subsequent waveform plotting

% %this has been done only the first time
% data_mat = zeros(204,length(timel),108,6);
% %%
% data_mat(:,:,:,deviant) = cat(3,dataetero,datahomo);
% %%
% save data_alldeviantsallsubjectseterohomo.mat data_mat

%% LBPD_startup_D

%starting up some of the functions that I wrote for the following analysis
pathl = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD'; %path to stored functions
addpath(pathl);
LBPD_startup_D(pathl);

%% reshaping data for later MCS calculation

%Please note that MCS involves stochastic processes, therefore you may see slightly different numbers than the ones reported in the Paper. However, these small differences do not affect the results in any way.
deviant = 6; %set the deviant you want to compute statistics for (1 = intensity; 2 = location; 3 = pitch; 4 = rhythm; 5 = slide; 6 = timbre; 7 = average over all deviants)
contrast = 1; %set 1 for cond1 > cond2; set 0 for cond2 > cond1 (this should be done for gradiometers only!)

%input data
pthresh = 0.01; %p-value threshold for data
%set time-points extremes (default should be 1 and 92)
mint = 1; %1 corresponds to time = 0 seconds, so the onset of the deviant stimulus (trials have been previously considered with baseline and baseline-corrected, then we extracted only the time-window of interest for MMN)
maxt = 92;
grad_lab = 1; %set 1 to grad only; 2 for mag only; 3 for both (we computed analysis for gradiometers only)

%actual computation
list = dir('/scratch5/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/COMT_MMN/MEGSensors/*mat');
%loading input data
load([list(1).folder '/' list(deviant+1).name]);
load('/scratch5/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/COMT_MMN/MEGSensors/ztime.mat') %loading time vector (time in seconds)
%binarizing p-values
P = DATAP; %trick for looking only into the positive t-values (this is for avoiding to mix significances coming from etero > homo and homo > etero at the same time)
P(P < pthresh) = 1; %binarizing p-values according to threshold
P(P < 1) = 0;
% TSTAT_mag = P(103:204,mint:maxt);
%gradiometers
TSTAT_grad = P(1:102,mint:maxt);
if contrast == 1
    TSTAT_grad(TSTAT(1:102,mint:maxt)<0) = 0; %taking only p-values corresponding to positive t-values (cond1 > cond2)
else
    TSTAT_grad(TSTAT(1:102,mint:maxt)>0) = 0; %taking only p-values corresponding to positive t-values (cond1 > cond2)
end
%positive magnetometers
TSTAT_mag_pos = P(103:204,mint:maxt);
TSTAT_mag_pos(TSTAT(103:204,mint:maxt)<0) = 0; %taking only p-values corresponding to positive t-values (cond1 > cond2)
%negative magnetometers
TSTAT_mag_neg = P(103:204,mint:maxt);
TSTAT_mag_neg(TSTAT(103:204,mint:maxt)>0) = 0; %taking only p-values corresponding to negative t-values (cond1 > cond2)
load('/scratch5/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/COMT_MMN/MEGSensors/zlab.mat')
%load a decent 2D approximation in a matrix of the MEG channels location
[~,~,raw_channels] = xlsread('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/MatrixMEGChannelLayout_With0_2.xlsx');
%actual input structure and function for reshaping data
S2 = [];
S2.label = label;
S2.TSTAT_mag_pos = TSTAT_mag_pos;
S2.TSTAT_mag_neg = TSTAT_mag_neg;
S2.TSTAT_grad = TSTAT_grad;
S2.raw_channels = raw_channels;
[MAG_data_pos, MAG_data_neg, GRAD_data] = MEG_sensors_MCS_reshapingdata_LBPD_D(S2);
%actual input structure and function for Monte Carlo simulations
S2 = [];
%gradiometers data
S2.data(:,:,:,1) = GRAD_data;
%magnetometers data
S2.data(:,:,:,2) = MAG_data_pos;
S2.data(:,:,:,3) = MAG_data_neg;
S2.sensortype = grad_lab;
S2.MEGlayout = cell2mat(raw_channels);
S2.permut = 1000;
S2.clustmax = 1;
S2.permthresh = 0.002;

[MAG_clust_pos, MAG_clust_neg, GRAD_clust] = MEG_sensors_MonteCarlosim_LBPD_D(S2); %actual MCS function

%% plotting significant clusters (topoplot)

for ii = 1:size(GRAD_clust,1)
    S2 = [];
    S2.PPall(1) = {GRAD_clust}; S2.PPall(2) = {{[]}}; S2.PPall(3) = {{[]}};
    S2.fieldtrip_mask = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External';
    S2.clustnum_grad = ii;
    S2.zlim = [0 21];
    S2.time = time(31:122);
    
    MEG_sensors_MCS_plottingclusters_LBPD_D(S2)
end
set(gcf,'Color','w')
%colormap with white for 0 valuescd /pro
x = []; x.bottom = [0 0 0.5]; x.botmiddle = [0 0.5 1]; x.middle = [1 1 1]; x.topmiddle = [1 0 0]; x.top = [0.6 0 0]; %red - blue
colormap(bluewhitered_PD(0,x))

%% plotting significant clusters (waveforms)

clustnum = 2; %cluster number that you want to plot (leave it to 2 if you want to have the cluster significant for all deviants together)
condition = 1:6; %plot single condition (e.g. S.condition = 3) or average waveform across conditions (e.g. S.condition = 2:7) (1 = intensity; 2 = location; 3 = pitch; 4 = rhythm; 5 = slide; 6 = timbre)

%additional inputs (that should not be changed) and actual computation
load('/scratch5/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/COMT_MMN/MEGSensors/zdata_alldeviantsallsubjectseterohomo.mat'); %loading data properly reshaped for the following function
load('/scratch5/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/COMT_MMN/MEGSensors/ztime.mat') %loading time vector (time in seconds)
load('/scratch5/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/COMT_MMN/MEGSensors/zchanlabels.mat'); %loading chanlabels coherently with the order of the data
load('/scratch5/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/COMT_MMN/MEGSensors/grad_clust_indices.mat'); %loading the channels forming the main significant cluster across all deviants

%structure for the function
S = [];
S.data = data_mat;
S.condition_n = condition;
S.groups = {'hetero','homo'}; %Group labels
gsubj{1} = [1:58]; gsubj{2} = [59:108];
S.gsubj = gsubj;
S.signtp = [0.020 0.057];
S.time_real = time;
S.legendl = 0;
S.colorline = {'r', 'b'}; %Select colorline for each group
S.sensors = 0; % -1 for magnetometer, 0 for gradiometers (gradiometers have already been combined)
S.x_lim = []; % Set x limits
S.y_lim = [1.0e-13 2.7e-12]; %Set y limits
S.conditions = {'Intensity','Localization','Pitch','Rhythm','Slide','Timbre'};
S.chanlabels = chanlabels;
S.STE = 1;
for kk = clustnum%1:size(GRAD_clust,1)
    S.data = data_mat; %data to be plotted
    clustplot = GRAD_clust{kk,3}; %slightly elaborated way to get the channel original IDs of channels forming the significant cluster that you are considering
    jes = [];
    for ii = 1:size(clustplot,1)
        jes(ii) = find(cellfun(@isempty,strfind(chanlabels,clustplot{ii,1})) == 0);
    end
    S.chans_index = jes; %plot waveform at single channel or the average across multiple channels (specify channel index - you can find the channel labels in 'chanlabels'); set to 0 to plot all channels
    plot_sensors_wavebis2(S) %actual function
end

%% plotting 3 groups according to COMT polymorphism

condition = 1:6; %plot single condition (e.g. S.condition = 3) or average waveform across conditions (e.g. S.condition = 2:7)
%(1 = intensity; 2 = location; 3 = pitch; 4 = rhythm; 5 = slide; 6 = timbre)
clustnum = 5; %cluster number that you want to plot (leave it to 2 if you want to have the cluster significant for all deviants together)

%additional inputs (that should not be changed) and actual computation
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/MumufeBDNF')
if length(condition) ~= 1
    load('/scratch5/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/COMT_MMN/MEGSensors/grad_clust_indices.mat'); %loading the channels forming the main significant cluster across all deviants
    clustnum = 2; %cluster number that you want to plot (leave it to 2 if you want to have the cluster significant for all deviants together)
end
load('/scratch5/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/COMT_MMN/MEGSensors/ztime.mat') %loading time vector (time in seconds)
load('/scratch5/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/COMT_MMN/MEGSensors/zchanlabels.mat'); %loading chanlabels coherently with the order of the data
load('/scratch5/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/COMT_MMN/MEGSensors/zdata_alldeviantsallsubjectseterohomo.mat'); %loading data properly reshaped for the following function
%structure for the function
S = [];
S.data = data_mat;
S.condition_n = condition;
% S.groups = {'hetero','homo'}; %Group labels
S.groups = {'hetero','met/met','val/val'}; %Group labels
% subjs_hetero = {'003','005','011','012','015','017','019','029','031','032','036','037','038','039','040','043','046','047','049','050','054','055','056','058','059','061','062','063','064','065','066','068','072','079','080','083','084','090','092','093','095','102','107','108','109','110','112','121','123','124','128','130','133','134','135','136','137','139'};
% subjs_homo = {'002','004','007','014','021','022','025','026','027','028','030','033','034','035','045','048','051','052','053','057','067','069','070','071','073','074','075','076','078','081','082','085','086','097','103','105','106','111','113','114','116','117','118','120','122','125','131','132','138','140'};
subjs_2met = {'002','004','007','014','022','025','026','028','045','048','053','057','067','069','070','071','075','076','081','086','097','103','111','113','114','116','118','120','125','131','132'}; %correspond to '0' in the .csv file'098' is missing
subjs_2val = {'021','027','030','033','034','035','051','052','073','074','078','082','085','105','106','117','122','138','140'}; %'2' in the .csv file; 13,42,73,74,119 missing
% gsubj{1} = [1:58]; gsubj{2} = [59:108];
gsubj{1} = [1:58]; gsubj{2} = [59:62 64:66 68 73:74 77:82 85:86 88 91:93 96:99 101:102 104:106]; gsubj{3} = [63 67 69:72 75:76 83:84 87 89:90 94:95 100 103 107:108];
S.gsubj = gsubj;
S.signtp = [];
S.time_real = time;
S.legendl = 0;
S.colorline = {'r','b',[0.4 0.4 0.4]}; %Select colorline for each group
S.sensors = 0; % -1 for magnetometer, 0 for gradiometers (gradiometers have already been combined)
S.x_lim = []; % Set x limits
S.y_lim = []; %Set y limits
S.conditions = {'Intensity','Localization','Pitch','Rhythm','Slide','Timbre'};
S.chanlabels = chanlabels;
S.STE = 1;
for kk = clustnum%1:size(GRAD_clust,1)
    S.data = data_mat; %data to be plotted
    clustplot = GRAD_clust{kk,3}; %slightly elaborated way to get the channel original IDs of channels forming the significant cluster that you are considering
    jes = [];
    for ii = 1:size(clustplot,1)
        jes(ii) = find(cellfun(@isempty,strfind(chanlabels,clustplot{ii,1})) == 0);
    end
    S.chans_index = jes; %plot waveform at single channel or the average across multiple channels (specify channel index - you can find the channel labels in 'chanlabels'); set to 0 to plot all channels
    plot_sensors_wavebis2(S) %actual function
end

%% trick to remove space from chanlabels.. useful for some computation but not needed now.. I left it here only to remember it
% 
% chanlabels2 = cell(1,204);
% for ii = 1:length(chanlabels)
%     idx = find(isspace(chanlabels{ii})==1); %removing space within data labels
%     chanlabels2{1,ii} = [chanlabels{ii}(1:idx-1) chanlabels{ii}(idx+1:end)]; %making data labels equal to the ones of topoplot function for later ICA usage..
% end
% chanlabels = chanlabels2;
% save chanlabels.mat chanlabels

%%

%% *** SOURCE RECONSTRUCTION ***

%% SOURCE RECONSTRUCTION (BEAMFORMING) AND STATISTICAL ANALYSIS

%% codes to get a "normal" oat for source reconstruction (so all subjects in the same folder)..

%settings for cluster (parallel computing)
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/scripts_osl_learningbach/Cluster')
clusterconfig('slot', 1); % set manually the job cluster slots
% between 1 and 12 (n x 8gb of ram)
% clusterconfig('scheduler', 'none'); % set automatically the long run queue
clusterconfig('scheduler', 'cluster'); % set automatically the long run queue
clusterconfig('long_running', 1); % set automatically the long run queue
spm_list = dir('/scratch5/MINDLAB2017_MEG-LearningBach/Portis/es*');

% v = [64 70 136 142 146 154 230];
oat = [];
for ii = 2:2:length(spm_list)
    processed_file = [spm_list(ii).folder '/' spm_list(ii).name];
    oat.source_recon.D_epoched(ii/2)         = {processed_file};
    oat.source_recon.sessions_to_do(ii/2) = {[str2double(spm_list(ii).name(13:16))]}; %sessions to do among the file_list (subject 39 excluded because of problems during data collection..)
    oat.source_recon.results_fnames(ii/2) = {['session' num2str(str2double(spm_list(ii).name(13:16))) '_recon']};
    D_epoched = spm_eeg_load(processed_file);
    D_epoched = D_epoched.montage('switch',1); %switch the montage to 1 in order to be safer (so we have the AFRICA denoised data)
    D_epoched.save();
end
% Beamform
pca_dim_1 = 50;
oat.source_recon.pca_dim = pca_dim_1(ones(length(oat.source_recon.D_epoched),1)); %this seems to be necessary. It is for having the same number of values than the number of input D objects
oat.source_recon.modalities = {'MEGMAG'; 'MEGPLANAR'};
oat.source_recon.conditions        = {'Standard','Pitch','Timbre','Localization','Intensity','Slide','Rhythm'};
oat.source_recon.gridstep          = 8; % in mm
oat.source_recon.time_range        = [-0.1 0.39]; % time range in secs
oat.source_recon.freq_range        = [0.1 40]; % frequency range in Hz
oat.source_recon.type              = 'Scalar';
oat.source_recon.method            = 'beamform';
oat.source_recon.normalise_method  = 'mean_eig';
oat.source_recon.forward_meg       = 'Single Shell';
oat.source_recon.report.do_source_variance_maps = 1;
oat.source_recon.dirname           = [spm_list(1).folder '/sourcetryindsubj600Hz/firstlevel'];
   
jobid = job2cluster(@cluster_beamforming,oat); %running with parallel

%% moving files computed independently for each subject into a common new folder (trying to go back to the "usual" way of running oat analysis pipeline in OSL)..

%%% TO BE RUN ONLY THE FIRST TIME %%%

list = dir('/scratch5/MINDLAB2017_MEG-LearningBach/Portis/sourcetryindsubj600Hz/0*');
bs = '/scratch5/MINDLAB2017_MEG-LearningBach/Portis/sourcetryindsubj600Hz/firstlevel.oat';
for ii = 1:length(list)
    %concat file .dat
    as = ['/scratch5/MINDLAB2017_MEG-LearningBach/Portis/sourcetryindsubj600Hz/' list(ii).name '/concatMfsession' num2str(str2double(list(ii).name(1:4))) '_spm_meeg.dat'];
    status = movefile(as,bs)
    %concat file .mat
    as = ['/scratch5/MINDLAB2017_MEG-LearningBach/Portis/sourcetryindsubj600Hz/' list(ii).name '/concatMfsession' num2str(str2double(list(ii).name(1:4))) '_spm_meeg.mat'];
    status = movefile(as,bs)
    %session file
    as = ['/scratch5/MINDLAB2017_MEG-LearningBach/Portis/sourcetryindsubj600Hz/' list(ii).name '/session' num2str(str2double(list(ii).name(1:4))) '_recon.mat'];
    status = movefile(as,bs)
end

%% settings for cluster (parallel computing)

addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/scripts_osl_learningbach/Cluster')
clusterconfig('slot', 1); % set manually the job cluster slots
% between 1 and 12 (n x 8gb of ram)
% clusterconfig('scheduler', 'none'); % set automatically the long run queue
clusterconfig('scheduler', 'cluster'); % set automatically the long run queue
clusterconfig('long_running', 1); % set automatically the long run queue

%% first level (each experimental block for each subject, independently)

%FIRST LEVEL
design_matrix_summary = {};
design_matrix_summary{1} = [1 0 0 0 0 0 0];design_matrix_summary{2} = [0 1 0 0 0 0 0]; design_matrix_summary{3}=[0 0 1 0 0 0 0];design_matrix_summary{4}=[0 0 0 1 0 0 0]; design_matrix_summary{5}=[0 0 0 0 1 0 0]; design_matrix_summary{6}=[0 0 0 0 0 1 0]; design_matrix_summary{7}=[0 0 0 0 0 0 1];
oat.first_level.design_matrix_summary = design_matrix_summary;
% contrasts to be calculated:
oat.first_level.contrast = {};
% contrast design matrix
oat.source_recon.conditions = {'Standard','Pitch','Timbre','Localization','Intensity','Slide','Rhythm'};
oat.first_level.contrast{1} = [1 0 0 0 0 0 0]'; %Standard
oat.first_level.contrast{2} = [0 1 0 0 0 0 0]'; %Pitch
oat.first_level.contrast{3} = [0 0 1 0 0 0 0]'; %Timbre
oat.first_level.contrast{4} = [0 0 0 1 0 0 0]'; %Localization
oat.first_level.contrast{5} = [0 0 0 0 1 0 0]'; %Intensity
oat.first_level.contrast{6} = [0 0 0 0 0 1 0]'; %Slide
oat.first_level.contrast{7} = [0 0 0 0 0 0 1]'; %Rhythm
oat.first_level.contrast{8} = [-1 1 0 0 0 0 0]'; %Pitch - Standard
oat.first_level.contrast{9} = [-1 0 1 0 0 0 0]'; %Timbre - Standard
oat.first_level.contrast{10} = [-1 0 0 1 0 0 0]'; %Localization - Standard
oat.first_level.contrast{11} = [-1 0 0 0 1 0 0]'; %Intensity - Standard
oat.first_level.contrast{12} = [-1 0 0 0 0 1 0]'; %Slide - Standard
oat.first_level.contrast{13} = [-1 0 0 0 0 0 1]'; %Rhythm - Standard
%contrast names
oat.first_level.contrast_name = {};
oat.first_level.contrast_name{1} = 'Standard';
oat.first_level.contrast_name{2} = 'Pitch';
oat.first_level.contrast_name{3} = 'Timbre';
oat.first_level.contrast_name{4} = 'Localization';
oat.first_level.contrast_name{5} = 'Intensity';
oat.first_level.contrast_name{6} = 'Slide';
oat.first_level.contrast_name{7} = 'Rhythm';
oat.first_level.contrast_name{8} = 'Pitch - Standard';
oat.first_level.contrast_name{9} = 'Timbre - Standard';
oat.first_level.contrast_name{10} = '%Localization - Standard';
oat.first_level.contrast_name{11} = 'Intensity - Standard';
oat.first_level.contrast_name{12} = 'Slide - Standard';
oat.first_level.contrast_name{13} = 'Rhythm - Standard';
%original contrasts, try1 and try2 seem to give the exact same results
oat.first_level.report.first_level_cons_to_do = [1 2 3 4 5 6 7 8 9 10 11 12 13]; %[1 2 3]; %better to do 3 2 1 in order to get the information for the peak value of the contrast old-new
oat.first_level.time_range = [-0.1 0.39];
oat.first_level.post_tf_downsample_factor = 1;
%slow negativity
oat.first_level.cope_type = 'coape'; %this should calculate absolute values of glm copes before comparing different conditions (otherwise, by default the absolute value is done on the copes after comparing different conditions); the absolute value here seems to be necessary because of the sign ambiguity issue introduced by source reconstruction
%N100
% oat.first_level.cope_type = 'none'; %this should calculate absolute values of glm copes before comparing different conditions (otherwise, by default the absolute value is done on the copes after comparing different conditions); the absolute value here seems to be necessary because of the sign ambiguity issue introduced by source reconstruction
% oat.first_level.name = ['wholebrain_first_level_BC_nocoape'];
oat.first_level.name = ['wholebrain_first_level_BC'];
oat.first_level.bc = ones(1,13);
%to add if the oat has not been automatically saved
list = dir('/scratch5/MINDLAB2017_MEG-LearningBach/Portis/sourcetryindsubj600Hz/0*');
for ii = 1:length(list)
    oat.first_level.results_fnames(ii) = {['session' num2str(str2double(list(ii).name(1:4))) '_wholebrain_first_level_BC.mat']};
%     oat.first_level.results_fnames(ii) = {['session' num2str(str2double(list(ii).name(1:4))) '_wholebrain_first_level_BC_nocoape.mat']};
end

%% running first level on parallel computing
for ii = 1:length(list)
    oat.first_level.sessions_to_do = [];
    oat.first_level.sessions_to_do = [ii]; %here it seems that the session indexes oat.source_recon.results_fnames{ii} is directly related to the sequential 
%     oat.first_level.sessions_to_do = [str2double(list(ii).name(1:4))];
    jobid = job2cluster(@cluster_beamfirstlevel,oat);
end

%% SUBJECT LEVEL

%in this case it does not do anything since I have only one experimental
%block for each participant (however, for computational reasons, I have to
%run it)
%this is needed to read the proper subjects..
for ii = 1:140
    oat.first_level.results_fnames(ii) = {['session' num2str(ii) '_wholebrain_first_level_BC.mat']};
%     oat.first_level.results_fnames(ii) = {['session' num2str(ii) '_wholebrain_first_level_BC_nocoape.mat']};
end
% this is if you have a perfect correspondance between sessions and subjects
oat.subject_level.session_index_list = cell(1,140); %sarebbe length(subjects)
oat.subject_level.name = 'MMN';
oat.subject_level.subjects_to_do = [];
oat.subject_level.subjects_to_do = [1:140]; %crucial to then get proper subject IDs
%to update the names of the results files of subject level
for ii = 1:140
    oat.subject_level.session_index_list{ii} = ii;
    oat.subject_level.results_fnames(ii) = {['subject' num2str(ii) '_wholebrain_first_level_BC_MMN.mat']};
end

%% settings for cluster (parallel computing)

addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/scripts_osl_learningbach/Cluster')
clusterconfig('slot', 4); % set manually the job cluster slots
% between 1 and 12 (n x 8gb of ram)
% clusterconfig('scheduler', 'none'); % set automatically the long run queue
clusterconfig('scheduler', 'cluster'); % set automatically the long run queue
clusterconfig('long_running', 1); % set automatically the long run queue

%% GROUP LEVEL

%list of subjects divided according to COMT gene
subjs_homo = {'002','004','007','014','021','022','025','026','027','028','030','033','034','035','045','048','051','052','053','057','067','069','070','071','073','074','075','076','078','081','082','085','086','097','103','105','106','111','113','114','116','117','118','120','122','125','131','132','138','140'};
subjs_hetero = {'003','005','011','012','015','017','019','029','031','032','036','037','038','039','040','043','046','047','049','050','054','055','056','058','059','061','062','063','064','065','066','068','072','079','080','083','084','090','092','093','095','102','107','108','109','110','112','121','123','124','128','130','133','134','135','136','137','139'};

%list of subjects with source reconstruction
list_source = dir(['/scratch5/MINDLAB2017_MEG-LearningBach/Portis/sourcetryindsubj600Hz/0*']);
%making sure that there is correspondance between subjects divided according to COMT gene and subjects that had source reconstrcted results
%hetero..
et_s = [];
for ii = 1:length(subjs_hetero)
    for jj = 1:length(list_source) 
        if str2double(subjs_hetero{ii}) == str2double(list_source(jj).name(1:4))
            et_s = cat(1,et_s,str2double(subjs_hetero{ii}));
            break
        end
    end    
end
%homo..
hm_s = [];
for ii = 1:length(subjs_homo)
    for jj = 1:length(list_source) 
        if str2double(subjs_homo{ii}) == str2double(list_source(jj).name(1:4))
            hm_s = cat(1,hm_s,str2double(subjs_homo{ii}));
            break
        end
    end    
end
oat.group_level = [];
oat.group_level.name = 'group_level_eterovshomo_BC'; %OBS!! REMEMBER TO UPDATE THE NAME!
oat.group_level.subjects_to_do = [];
oat.group_level.subjects_to_do = cat(1,et_s,hm_s)';
%other usual settings
%results name
% oat.group_level.results_fnames = ['wholebrain_first_level_BC_MMN' '_' oat.group_level.name '.mat'];
oat.group_level.results_fnames = ['wholebrain_first_level_BC_MMN' '_' oat.group_level.name '.mat'];
% Spatial and temporal averaging options
oat.group_level.time_range = [-0.1 0.39];
oat.group_level.space_average = 0;
oat.group_level.time_average = 0;
oat.group_level.time_smooth_std = 0; % secs
oat.group_level.use_tstat = 0;
%path to AAL template
% oat.group_level.mask_fname = '/projects/MINDLAB2017_MEG-LearningBach/scripts/osl/AAL_Nifti/aal_8mm_try5.nii.gz';
% Spatial and temporal smoothing options
oat.group_level.spatial_smooth_fwhm = 0; % mm
oat.group_level.group_varcope_time_smooth_std = 100;
oat.group_level.group_varcope_spatial_smooth_fwhm = 100; % smooths the variance of the group copes. It is recommended to do this.
%store copes (useful for doing the permutation test later, otherwise it needs to compute again the group level analysis)
oat.group_level.store_lower_level_copes = 1;
% Set up design matrix and contrasts
%this is if you have only the general mean across the all participants
% oat.group_level.group_design_matrix =
% ones(1,length(oat.group_level.subjects_to_do)); %if you want all of the
% participantsoa
oat.group_level.group_contrast = [];
% oat.group_level.group_contrast{1} = [1];
% oat.group_level.group_contrast_name = {};
% oat.group_level.group_contrast_name{1} = 'mean';
oat.group_level.glm_method='fixed_effects'; %ols or fixed-effects
% Define which contrasts to perform for the report
oat.group_level.first_level_contrasts_to_do = [1,2,3,4,5,6,7,8,9,10,11,12,13]; % list of first level contrasts to run the group analysis on
oat.group_level.report.first_level_cons_to_do = [1]; % purely used to determine which contrasts are shown in the report, the first one in list determines the contrast used to find the vox, time, freq with the max stat
oat.group_level.report.group_level_cons_to_do = [1]; % purely used to determine which contrasts are shown in the report, the first one in list determines the contrast used to find the vox, time, freq with the max stat
oat.group_level.report.show_lower_level_copes = 0;
oat.group_level.report.show_lower_level_cope_maps = 0;
%with two groups
%1s for mus1 and 0s for mus0 for the design matrix
a2 = zeros(2,length(oat.group_level.subjects_to_do));
a2(1,1:58) = 1; %etero
a2(2,59:end) = 1; %homo
oat.group_level.group_design_matrix = a2;
%contrast specification
oat.group_level.group_contrast{1} = [1 -1]'; %group1 > group2, remember the '!!
oat.group_level.group_contrast_name{1} = 'Etero > Homo';
oat.group_level.group_contrast{2} = [-1 1]'; %group1 > group2, remember the '!!
oat.group_level.group_contrast_name{2} = 'Homo > Etero';

jobid = job2cluster(@cluster_beamgrouplevel,oat);

%% creatig nifti images with statistics..

% load('/scratch5/MINDLAB2017_MEG-LearningBach/Portis/sourcetryindsubj600Hz/firstlevel.oat/oat_wholebrain_first_level_BC_MMN_group_level_everybody_BC.mat');
S2 = [];
S2.oat = oat;
S2.stats_fname = oat.group_level.results_fnames;
S2.first_level_contrasts = 1:13; %remember that in this way you define the order of the output (this numbers refers to the order (numbers) defined in the contrasts; for example here tstat3 refers to the contrast old-new)
S2.group_level_contrasts = [1 2];
S2.resamp_gridstep = oat.source_recon.gridstep;

jobid = job2cluster(@cluster_oat_save_nii_stats,S2);

%% settings for cluster (parallel computing)

addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/scripts_osl_learningbach/Cluster')
clusterconfig('slot', 4); % set manually the job cluster slots
% between 1 and 12 (n x 8gb of ram)
% clusterconfig('scheduler', 'none'); % set automatically the long run queue
clusterconfig('scheduler', 'cluster'); % set automatically the long run queue
clusterconfig('long_running', 1); % set automatically the long run queue

%% cluster-based permutation test

load('/scratch5/MINDLAB2017_MEG-LearningBach/Portis/sourcetryindsubj600Hz/firstlevel.oat/oat_wholebrain_first_level_BC_MMN_group_level_eterovshomo_BC.mat')

S = [];
S.oat = oat;
% S.cluster_stats_thresh = 0.7;
% S.cluster_stats_thresh = 3.0; %tstat corresponding to p-value = .0017 (approximally .002) with 107 degrees of freedom
S.cluster_stats_thresh = 1.7; %tstat corresponding to p-value = .047
S.cluster_stats_nperms = 5000; % we normally recommend doing 5000 perms
S.first_level_copes_to_do = [9];
S.group_level_copes_to_do = [1,2];
S.group_varcope_spatial_smooth_fwhm = S.oat.group_level.group_varcope_spatial_smooth_fwhm;
S.write_cluster_script = 0;
% S.time_range = [0.153 0.280]; %slide vs standard
% S.time_range = [0.003 0.080]; %intensity vs standard
% S.time_range = [0.033 0.077]; %localization vs standard
% S.time_range = [0.263 0.306]; %pitch vs standard
% S.time_range = [0.226 0.276]; %rhythm vs standard
S.time_range = [0.020 0.057]; %timbre vs standard
S.time_average = 1;
% Run the permutations (on cluster.. parallel computing)
jobid = job2cluster(@clusterbasedpermutation_osl,S);

%[ gstats ] = oat_cluster_permutation_testing(S);

%The outputted results were then used as mask for plotting only the
%significant voxels and their corresponding statistics. This operation was
%done by using fslmath.

%We used the following example line in the terminal:
% fslmaths tstat3_gc1_8mm.nii.gz.nii.gz -mas clustere_corrp_tstat3_gc1_8mm.nii.gz contr.nii.gz

%After that, we used Workbench for producing the
%images that we presented in the paper.
%In the following section we report the images that we obtained and used to
%produce a detailed statistical output concerning the significant voxels.
%Additional clarification about these algorithms can be found in the paper.
%Furthermore, you are very welcome to contatct us if you have specific
%questions.

%% Runnging again source reconstruction (and statistics on it)

%This was done since we did not compute all contrasts in the first place..
%(we computed only contrasts for the 6 deviants independently and not for
%them together).
%This trivial oversight forced us to compute again the steps already reported above.
%To increase clarity, we decided to report in this script all steps again, as follows, so
%that readers can more easily see the pipeline employed.

%% codes to get a "normal" oat for source reconstruction (so all subjects in the same folder)..

spm_list = dir('/scratch5/MINDLAB2017_MEG-LearningBach/Portis/es*');
% v = [64 70 136 142 146 154 230];
oat = [];
for ii = 2:2:length(spm_list)
    processed_file = [spm_list(ii).folder '/' spm_list(ii).name];
    oat.source_recon.D_epoched(ii/2)         = {processed_file};
    oat.source_recon.sessions_to_do(ii/2) = {[str2double(spm_list(ii).name(13:16))]}; %sessions to do among the file_list (subject 39 excluded because of problems during data collection..)
    oat.source_recon.results_fnames(ii/2) = {['session' num2str(str2double(spm_list(ii).name(13:16))) '_recon']};
    D_epoched = spm_eeg_load(processed_file);
    D_epoched = D_epoched.montage('switch',1); %switch the montage to 1 in order to be safer (so we have the AFRICA denoised data)
    D_epoched.save();
end
% Beamform
pca_dim_1 = 50;
oat.source_recon.pca_dim = pca_dim_1(ones(length(oat.source_recon.D_epoched),1)); %this seems to be necessary. It is for having the same number of values than the number of input D objects
oat.source_recon.modalities = {'MEGMAG'; 'MEGPLANAR'};
oat.source_recon.conditions        = {'Standard','Pitch','Timbre','Localization','Intensity','Slide','Rhythm'};
oat.source_recon.gridstep          = 8; % in mm
oat.source_recon.time_range        = [-0.1 0.39]; % time range in secs
oat.source_recon.freq_range        = [0.1 40]; % frequency range in Hz
%     oat.source_recon.freq_range        = [0.1 10]; % frequency range in Hz
%S.source_recon.pca_order         = 250;
oat.source_recon.type              = 'Scalar';
oat.source_recon.method            = 'beamform';
oat.source_recon.normalise_method  = 'mean_eig';
oat.source_recon.forward_meg       = 'Single Shell';
%S.source_recon.prefix            = '';
oat.source_recon.report.do_source_variance_maps = 1;
oat.source_recon.dirname           = [spm_list(1).folder '/sourcetryindsubj600Hz/firstlevel'];

jobid = job2cluster(@cluster_beamforming,oat); %running with parallel

%% settings for cluster (parallel computing)

addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/scripts_osl_learningbach/Cluster')
clusterconfig('slot', 4); % set manually the job cluster slots
% between 1 and 12 (n x 8gb of ram)
% clusterconfig('scheduler', 'none'); % set automatically the long run queue
clusterconfig('scheduler', 'cluster'); % set automatically the long run queue
clusterconfig('long_running', 1); % set automatically the long run queue

%% first level (each experimental block for each subject, independently)

%FIRST LEVEL
design_matrix_summary = {};
design_matrix_summary{1} = [1 0 0 0 0 0 0];design_matrix_summary{2} = [0 1 0 0 0 0 0]; design_matrix_summary{3}=[0 0 1 0 0 0 0];design_matrix_summary{4}=[0 0 0 1 0 0 0]; design_matrix_summary{5}=[0 0 0 0 1 0 0]; design_matrix_summary{6}=[0 0 0 0 0 1 0]; design_matrix_summary{7}=[0 0 0 0 0 0 1];
oat.first_level.design_matrix_summary = design_matrix_summary;
% contrasts to be calculated:
oat.first_level.contrast = {};
% contrast design matrix
oat.source_recon.conditions        = {'Standard','Pitch','Timbre','Localization','Intensity','Slide','Rhythm'};
oat.first_level.contrast{1} = [-6 1 1 1 1 1 1]'; %All deviants vs standard
%contrast names
oat.first_level.contrast_name = {};
oat.first_level.contrast_name{1} = 'Deviants - Standard';
%original contrasts, try1 and try2 seem to give the exact same results
oat.first_level.report.first_level_cons_to_do = [1]; %[1 2 3]; %better to do 3 2 1 in order to get the information for the peak value of the contrast old-new
oat.first_level.time_range = [-0.1 0.39];
oat.first_level.post_tf_downsample_factor = 1;
%slow negativity
oat.first_level.cope_type = 'coape'; %this should calculate absolute values of glm copes before comparing different conditions (otherwise, by default the absolute value is done on the copes after comparing different conditions); the absolute value here seems to be necessary because of the sign ambiguity issue introduced by source reconstruction
%N100
% oat.first_level.cope_type = 'none'; %this should calculate absolute values of glm copes before comparing different conditions (otherwise, by default the absolute value is done on the copes after comparing different conditions); the absolute value here seems to be necessary because of the sign ambiguity issue introduced by source reconstruction
oat.first_level.name = ['wholebrain_first_level_BC_alldeviantsminusstrandard']; %REMEMBER TO CHECK THIS NAME!!
oat.first_level.bc = ones(1,1);
%to add if the oat has not been automatically saved
list = dir('/scratch5/MINDLAB2017_MEG-LearningBach/Portis/sourcetryindsubj600Hz/0*');
for ii = 1:length(list)
    oat.first_level.results_fnames(ii) = {['session' num2str(str2double(list(ii).name(1:4))) '_wholebrain_first_level_BC_alldeviantsminusstrandard.mat']};
%     oat.first_level.results_fnames(ii) = {['session' num2str(str2double(list(ii).name(1:4))) '_wholebrain_first_level_BC_nocoape.mat']};
end

%% running first level on parallel computing
for ii = 1:length(list)
    oat.first_level.sessions_to_do = [];
    oat.first_level.sessions_to_do = [ii]; %here it seems that the session indexes oat.source_recon.results_fnames{ii} is directly related to the sequential 
%     oat.first_level.sessions_to_do = [str2double(list(ii).name(1:4))];
    jobid = job2cluster(@cluster_beamfirstlevel,oat);
end

%% SUBJECT LEVEL

%in this case it does not do anything since I have only one experimental
%block for each participant (however, for computational reasons, I have to
%run it)
%this is needed to read the proper subjects..
for ii = 1:140
    oat.first_level.results_fnames(ii) = {['session' num2str(ii) '_wholebrain_first_level_BC_alldeviantsminusstrandard.mat']};
%     oat.first_level.results_fnames(ii) = {['session' num2str(ii) '_wholebrain_first_level_BC_nocoape.mat']};
end
% this is if you have a perfect correspondance between sessions and subjects
oat.subject_level.session_index_list = cell(1,140); %sarebbe length(subjects)
oat.subject_level.name = 'MMN';
oat.subject_level.subjects_to_do = [];
oat.subject_level.subjects_to_do = [1:140]; %crucial to then get proper subject IDs
%to update the names of the results files of subject level
for ii = 1:140
    oat.subject_level.session_index_list{ii} = ii;
    oat.subject_level.results_fnames(ii) = {['subject' num2str(ii) '_wholebrain_first_level_BC_alldeviantsminusstrandard_MMN.mat']};
end

oat = osl_check_oat(oat);
oat.to_do = [0 0 1 0];
oat = osl_run_oat(oat);

%% settings for cluster (parallel computing)

addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/scripts_osl_learningbach/Cluster')
clusterconfig('slot', 4); % set manually the job cluster slots
% between 1 and 12 (n x 8gb of ram)
% clusterconfig('scheduler', 'none'); % set automatically the long run queue
clusterconfig('scheduler', 'cluster'); % set automatically the long run queue
clusterconfig('long_running', 1); % set automatically the long run queue


%% GROUP LEVEL

%list of subjects divided according to COMT gene
subjs_homo = {'002','004','007','014','021','022','025','026','027','028','030','033','034','035','045','048','051','052','053','057','067','069','070','071','073','074','075','076','078','081','082','085','086','097','103','105','106','111','113','114','116','117','118','120','122','125','131','132','138','140'};
subjs_hetero = {'003','005','011','012','015','017','019','029','031','032','036','037','038','039','040','043','046','047','049','050','054','055','056','058','059','061','062','063','064','065','066','068','072','079','080','083','084','090','092','093','095','102','107','108','109','110','112','121','123','124','128','130','133','134','135','136','137','139'};

%list of subjects with source reconstruction
list_source = dir(['/scratch5/MINDLAB2017_MEG-LearningBach/Portis/sourcetryindsubj600Hz/0*']);
%making sure that there is correspondance between subjects divided according to COMT gene and subjects that had source reconstrcted results
%hetero..
et_s = [];
for ii = 1:length(subjs_hetero)
    for jj = 1:length(list_source) 
        if str2double(subjs_hetero{ii}) == str2double(list_source(jj).name(1:4))
            et_s = cat(1,et_s,str2double(subjs_hetero{ii}));
            break
        end
    end    
end
%homo..
hm_s = [];
for ii = 1:length(subjs_homo)
    for jj = 1:length(list_source) 
        if str2double(subjs_homo{ii}) == str2double(list_source(jj).name(1:4))
            hm_s = cat(1,hm_s,str2double(subjs_homo{ii}));
            break
        end
    end    
end
oat.group_level = [];
oat.group_level.name = 'group_level_eterovshomo_BC2'; %OBS!! REMEMBER TO UPDATE THE NAME!
oat.group_level.subjects_to_do = [];
oat.group_level.subjects_to_do = cat(1,et_s,hm_s)';
%other usual settings
%results name
% oat.group_level.results_fnames = ['wholebrain_first_level_BC_MMN' '_' oat.group_level.name '.mat'];
oat.group_level.results_fnames = ['wholebrain_first_level_BC_alldeviantsminusstrandard_MMN' '_' oat.group_level.name '.mat'];
% Spatial and temporal averaging options
oat.group_level.time_range = [-0.1 0.39];
oat.group_level.space_average = 0;
oat.group_level.time_average = 0;
oat.group_level.time_smooth_std = 0; % secs
oat.group_level.use_tstat = 0;
%path to AAL template
% oat.group_level.mask_fname = '/projects/MINDLAB2017_MEG-LearningBach/scripts/osl/AAL_Nifti/aal_8mm_try5.nii.gz';
% Spatial and temporal smoothing options
oat.group_level.spatial_smooth_fwhm = 0; % mm
oat.group_level.group_varcope_time_smooth_std = 100;
oat.group_level.group_varcope_spatial_smooth_fwhm = 100; % smooths the variance of the group copes. It is recommended to do this.
%store copes (useful for doing the permutation test later, otherwise it needs to compute again the group level analysis)
oat.group_level.store_lower_level_copes = 1;
% Set up design matrix and contrasts
%this is if you have only the general mean across the all participants
% oat.group_level.group_design_matrix =
% ones(1,length(oat.group_level.subjects_to_do)); %if you want all of the
% participantsoa
oat.group_level.group_contrast = [];
% oat.group_level.group_contrast{1} = [1];
% oat.group_level.group_contrast_name = {};
% oat.group_level.group_contrast_name{1} = 'mean';
oat.group_level.glm_method='fixed_effects'; %ols or fixed-effects
% Define which contrasts to perform for the report
oat.group_level.first_level_contrasts_to_do = [1]; % list of first level contrasts to run the group analysis on
oat.group_level.report.first_level_cons_to_do = [1]; % purely used to determine which contrasts are shown in the report, the first one in list determines the contrast used to find the vox, time, freq with the max stat
oat.group_level.report.group_level_cons_to_do = [1]; % purely used to determine which contrasts are shown in the report, the first one in list determines the contrast used to find the vox, time, freq with the max stat
oat.group_level.report.show_lower_level_copes = 0;
oat.group_level.report.show_lower_level_cope_maps = 0;
%with two groups
%1s for mus1 and 0s for mus0 for the design matrix
a2 = zeros(2,length(oat.group_level.subjects_to_do));
a2(1,1:58) = 1; %etero
a2(2,59:end) = 1; %homo
oat.group_level.group_design_matrix = a2;
%contrast specification
oat.group_level.group_contrast{1} = [1 -1]'; %group1 > group2, remember the '!!
oat.group_level.group_contrast_name{1} = 'Etero > Homo';
oat.group_level.group_contrast{2} = [-1 1]'; %group1 > group2, remember the '!!
oat.group_level.group_contrast_name{2} = 'Homo > Etero';

jobid = job2cluster(@cluster_beamgrouplevel,oat);

%% cluster-based permutation test

load('/scratch5/MINDLAB2017_MEG-LearningBach/Portis/sourcetryindsubj600Hz/firstlevel.oat/oat_wholebrain_first_level_BC_alldeviantsminusstrandard_MMN_group_level_eterovshomo_BC2.mat')

S = [];
S.oat = oat;
% S.cluster_stats_thresh = 0.7;
% S.cluster_stats_thresh = 3.0; %tstat corresponding to p-value = .0017 (approximally .002) with 107 degrees of freedom
S.cluster_stats_thresh = 1.7; %tstat corresponding to p-value = .047
S.cluster_stats_nperms = 5000; % we normally recommend doing 5000 perms
S.first_level_copes_to_do = [1];
S.group_level_copes_to_do = [1,2];
S.group_varcope_spatial_smooth_fwhm = S.oat.group_level.group_varcope_spatial_smooth_fwhm;
S.write_cluster_script = 0;
S.time_range = [0.19 0.29]; %all deviants vs standard
S.time_average = 1;
% Run the permutations (on cluster.. parallel computing)
jobid = job2cluster(@clusterbasedpermutation_osl,S);
%[ gstats ] = oat_cluster_permutation_testing(S);

%The outputted results were then used as mask for plotting only the
%significant voxels and their corresponding statistics. This operation was
%done by using fslmath.

%We used the following example line in the terminal:
% fslmaths tstat3_gc1_8mm.nii.gz.nii.gz -mas clustere_corrp_tstat3_gc1_8mm.nii.gz contr.nii.gz

%After that, we used Workbench for producing the
%images that we presented in the paper.
%In the following section we report the images that we obtained and used to
%produce a detailed statistical output concerning the significant voxels.
%Additional clarification about these algorithms can be found in the paper.
%Furthermore, you are very welcome to contatct us if you have specific
%questions.


%%


%% creating a file with information on significant clusters/voxels outputted by the permutation test (THINK TO MAKE THIS A FUNCTION)

deviant = 1; %select the deviant that you want (1 = intensity; 2 = location; 3 = pitch; 4 = rhythm; 5 = slide; 6 = timbre; 7 = statistics for the average over the six deviants)

%path to file nifti
pathnii = dir('/scratch5/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/COMT_MMN/MEGSources/*gz'); %THIS MUST BE YOUR OWN DIRECTORY
%actual name plus path
fname = [pathnii(deviant).folder '/' pathnii(deviant).name];
%getting MNI coordinates of significant voxels within the provided image
[ mni_coords, xform ] = osl_mnimask2mnicoords(fname);
%loading the image
V = nii.load(fname);
%extracting statistics
VV = V(V~=0);
%indices of non-zero values of nifti image
VI = find(V~=0);
%path to AAL template
parcelfile = '/projects/MINDLAB2017_MEG-LearningBach/scripts/osl/AAL_Nifti/aal_8mm_try5.nii.gz';
%loading AAL labels
load('/projects/MINDLAB2017_MEG-LearningBach/scripts/scripts_osl_learningbach/LEiDA/AAL_labels.mat');
%extracting AAL coordinates information
K = nii.load(parcelfile);
%finding AAL non-zero coordinates values
% KI = find(K~=0);
%sorting results in order to have strongest voxels at the top (positive
%t-values) or at the bottom (negative t-values)
[VV2, II] = sort(VV,'descend');
VI = VI(II);
mni_coords = mni_coords(II,:);
%final cell
PD = cell(length(VV2),4);
%getting AAL indices
ROI = zeros(length(VI),1);
cnt = 0;
for ii = 1:length(VI)
    ROI(ii) = K(VI(ii));
    if ROI(ii) > 0 && ROI(ii) < 91
        cnt = cnt + 1;
        PD(cnt,1) = {lab(ROI(ii),3:end)}; %storing ROI
        PD(cnt,4) = {mni_coords(ii,:)}; %storing MNI coordinates
        if mni_coords(ii,1) > 0 %storing hemisphere
            PD(cnt,2) = {'R'};
        else
            PD(cnt,2) = {'L'};
        end
        PD(cnt,3) = {round(VV2(ii),2)}; %storing t-statistics
    end
end
PDn = cell2table(PD(~any(cellfun('isempty',PD),2),:)); %remove the possible empty cell
writetable(PDn,[pathnii(deviant).name(1:8) '.xlsx'],'Sheet',1)

%%

%% REVIEWERS' REQUEST -- PLOTTING STANDARD AND DEVIANT RESPONSES INDEPENDENTLY AND SEE WHETHER THEY ARE DIFFERENT IN RELATION TO COMT (INSTEAD OF HAVING THE PROPER MMN CALCULATED)

clustnum = 2; %cluster number that you want to plot (leave it to 2 if you want to have the cluster significant for all deviants together)
condition = 1; %plot single condition (e.g. S.condition = 3) or average waveform across conditions (e.g. S.condition = 2:7) (1 = intensity; 2 = location; 3 = pitch; 4 = rhythm; 5 = slide; 6 = timbre)
S.conditions = {'Intensity','Localization','Pitch','Rhythm','Slide','Timbre'};

%additional inputs (that should not be changed) and actual computation
load('/scratch5/MINDLAB2017_MEG-LearningBach/Portis/MEGsensor/sensor_data.mat'); %loading data properly reshaped for the following function
% load('/scratch5/MINDLAB2017_MEG-LearningBach/Portis/MEGsensor/MMNall_standardsubtracted.mat'); %loading data properly reshaped for the following function
load('/scratch5/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/COMT_MMN/MEGSensors/ztime.mat') %loading time vector (time in seconds)
% load('/scratch5/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/COMT_MMN/MEGSensors/zchanlabels.mat'); %loading chanlabels coherently with the order of the data
load('/scratch5/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/COMT_MMN/MEGSensors/grad_clust_indices.mat'); %loading the channels forming the main significant cluster across all deviants
%structure for the function
S = [];
S.data = data_mat;
S.condition_n = condition;
S.groups = {'hetero','homo'}; %Group labels
subjs_hetero = {'003','005','011','012','015','017','019','029','031','032','036','037','038','039','040','043','046','047','049','050','054','055','056','058','059','061','062','063','064','065','066','068','072','079','080','083','084','090','092','093','095','102','107','108','109','110','112','121','123','124','128','130','133','134','135','136','137','139'};
subjs_homo = {'002','004','007','014','021','022','025','026','027','028','030','033','034','035','045','048','051','052','053','057','067','069','070','071','073','074','075','076','078','081','082','085','086','097','103','105','106','111','113','114','116','117','118','120','122','125','131','132','138','140'};
listmeg = dir('/scratch5/MINDLAB2017_MEG-LearningBach/Portis/sourcetryindsubj600Hz/0*');
cnte = 0;
cnth = 0;
for ii = 1:length(listmeg)
    for pp = 1:length(subjs_hetero)
        if strcmp(listmeg(ii).name(2:4),subjs_hetero{pp})
            cnte = cnte + 1;
            ghet(cnte) = ii;
        end
    end
    for pp = 1:length(subjs_homo)
        if strcmp(listmeg(ii).name(2:4),subjs_homo{pp})
            cnth = cnth + 1;
            ghom(cnth) = ii;
        end
    end
end
gsubj{1} = ghet;
gsubj{2} = ghom;
S.gsubj = gsubj;
S.signtp = [];
S.time_real = time_sel;
S.legendl = 1;
S.colorline = {'r', 'b'}; %Select colorline for each group
S.sensors = 0; % -1 for magnetometer, 0 for gradiometers (gradiometers have already been combined)
S.x_lim = []; % Set x limits
S.y_lim = []; %Set y limits
S.conditions = {'Standard','Intensity','Localization','Pitch','Rhythm','Slide','Timbre'};
S.chanlabels = chanlabels;
S.STE = 2;
for kk = clustnum%1:size(GRAD_clust,1)
    S.data = data_mat; %data to be plotted
    clustplot = GRAD_clust{kk,3}; %slightly elaborated way to get the channel original IDs of channels forming the significant cluster that you are considering
    jes = [];
    for ii = 1:size(clustplot,1)
        jes(ii) = find(cellfun(@isempty,strfind(chanlabels,clustplot{ii,1})) == 0);
    end
    S.chans_index = jes; %plot waveform at single channel or the average across multiple channels (specify channel index - you can find the channel labels in 'chanlabels'); set to 0 to plot all channels
    plot_sensors_wavebis2(S) %actual function
end

%%