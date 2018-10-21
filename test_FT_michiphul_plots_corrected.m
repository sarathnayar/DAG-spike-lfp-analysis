% function test_FT_spike_field_coherence
% http://www.fieldtriptoolbox.org/tutorial/spikefield

cfg = [];
cfg.numtrl	= 20;
cfg.fsample     = 1000; % Hz
cfg.trllen      = 1; % s
        
cfg.s1.freq     = 10; % Hz
cfg.s1.phase    = 'random'; % 'random' or 0
cfg.s1.ampl     = 1;

cfg.s2.freq = 25; % Hz
cfg.s2.ampl = 0.5;
cfg.s2.phase = 0;

cfg.s3.freq = 60; % Hz
cfg.s3.ampl = 0.2;
cfg.s3.phase = 0;

cfg.noise.ampl = 0.2;

cfg.method	= 'superimposed';
cfg.output	= 'all'; % mixed or all

lfpcfg = cfg;


% In the method 'superimposed' the signal contains just the sum of the different frequency contributions:
%     s1: first frequency
%     s2: second frequency
%     s3: third frequency
% and the output consists of the following channels:
%     1st channel: mixed signal = s1 + s2 + s3 + noise
%     2nd channel: s1
%     3rd channel: s2
%     4th channel: s3
%     5th channel: noise
   
lfp_freq = [cfg.s1.freq cfg.s2.freq cfg.s3.freq];
lfp_amp =  [cfg.s1.ampl cfg.s2.ampl cfg.s3.ampl];


% spikes
spikeRates = 5:5:100; % Hz, overall spike rate
spikePhaseFreq = 1; % 1 or 2 or 3, s1 or s2 or s3 components
spikePhaseMean = [-pi]; % align to trough
spikePhaseStd  = [pi/60]; % rad, spread around spikePhaseMean
spikeLockProbs = 0:0.1:1; % probability of spike to be locked, the phase-locked firing rate would be max(spikeRate*spikeLockProb,lfp_freq(spikePhaseFreq))

% initialize variables to store ppc and coh
ppc = zeros(length(spikeLockProbs), length(spikeRates));
coh = zeros(length(spikeLockProbs), length(spikeRates));
% loop through spike locking probabilities
for i = 1:numel(spikeLockProbs)
    spikeLockProb = spikeLockProbs(i);
    disp(sprintf('Spike locking probability = %0.2f', spikeLockProb));
    % loop through different spike rates
    for j = 1:numel(spikeRates)
        spikeRate = spikeRates(j);
        disp(sprintf('Spike rate = %0.2f', spikeRate));
        cfg.numtrl	= 20;
        cfg.fsample     = 1000; % Hz
        cfg.trllen      = 1; % s
        data = ft_freqsimulation(lfpcfg);

        lockedSpikes = zeros(cfg.numtrl,cfg.trllen*cfg.fsample);
        Spikes = test_FT_simulate_spike_train(spikeRate*(1-spikeLockProb),cfg.trllen,cfg.numtrl); % non-locked spikes
        for t = 1:cfg.numtrl,
            % for each trial, find the phase (of the trough) of the corresponding freq component s
             s = data.trial{t}(spikePhaseFreq+1,:);
             taxis = data.time{t};
             p = angle(hilbert(s));
             trough_idx = find(diff(p)<-6); % find indices of the troughs

             lockedSpikes_idx = trough_idx(rand(size(trough_idx)) <= spikeRate*spikeLockProb/lfp_freq(spikePhaseFreq));

             % add some variability to relative phase of spikes and troughs
             lockedSpikes_idx = fix(lockedSpikes_idx+randn(size(lockedSpikes_idx))*spikePhaseStd*(cfg.fsample/lfp_freq(spikePhaseFreq)));
             lockedSpikes_idx = ig_limit_range_min_max(lockedSpikes_idx,1,cfg.fsample*cfg.trllen);
             lockedSpikes(t,lockedSpikes_idx)=1;
             Spikes(t,lockedSpikes_idx) = 1; % add locked spikes to spikes
             data.trial{t} = [data.trial{t}; Spikes(t,:)];

        end

        disp(sprintf('mean rate %.2f, mean locked rate %.2f %d spikes',spikeRate,spikeLockProb,mean(sum(Spikes,2)/cfg.trllen),mean(sum(lockedSpikes,2)/cfg.trllen),sum(sum(Spikes))));


        data.label = {'lfp1', 's1', 's2', 's3', 'noise', 'spk1'};
        % data_all = ft_appendspike([],data, Spikes)

        % Spike-triggered spectrum

        % method 1
        cfg           = [];
        cfg.method    = 'mtmconvol';
        % set the frequency to be the frequency of first LFP component
        cfg.foi       = 5:5:20;
        cfg.t_ftimwin = 5./cfg.foi; % 5 cycles per frequency
        cfg.taper     = 'hanning';
        cfg.spikechannel = data.label{6};
        cfg.channel      = data.label{1};
        stsConvol     = ft_spiketriggeredspectrum(cfg, data);

        % compute the statistics on the phases
        cfg               = [];
        cfg.method        = 'ppc0'; % compute the Pairwise Phase Consistency, can also be plv, ang, ppc1, ppc2, ral
        cfg.spikechannel  = stsConvol.label{1};
        cfg.channel       = stsConvol.lfplabel(1); % selected LFP channels
        cfg.avgoverchan   = 'unweighted'; % weight spike-LFP phases irrespective of LFP power
        cfg.timwin        = 'all'; % compute over all available spikes in the window
        cfg.latency       = 'maxperiod'; % [0 nanmax(stsConvol.trialtime(:))]; %
        statSts           = ft_spiketriggeredspectrum_stat(cfg,stsConvol);

        % get the ppc at frequency of first LFP component, 10 Hz
        ppc(i, j) = statSts.ppc0(statSts.freq == lfp_freq(1));


        % now assess coherence using Cronux
        
        params.Fs	=1000; % sampling frequency
        params.fpass	=[1 lfp_freq(1)]; % band of frequencies to be kept
        params.tapers	=[2 4]; % taper parameters
        params.pad	=2; % pad factor for fft
        params.err	=[2 0.05];
        params.trialave =1; % average coherence over trials

        % LFP: time x trials
        cfg			= [];
        cfg.trials		= 'all';
        cfg.channel		= data.label{1};
        data_lfp		= ft_selectdata(cfg, data);
        chr_data1		= test_FT_fieldtrip2chronux(data_lfp,'lfp');

        cfg.channel		= data.label{6};
        data_spikes		= ft_selectdata(cfg, data);
        chr_data2		= test_FT_fieldtrip2chronux(data_spikes,'spikes');

        [C,phi,S12,S1,S2,f,zerosp,confC,phistd,Cerr]=coherencycpt(chr_data1,chr_data2,params,0);
        
        % get the coh value at frequency of first LFP component
        coh(i, j) = C(end);
        
    end
    
end

% plots
rand_rgb = rand(numel(spikeLockProbs),3);
for i = 1:numel(spikeLockProbs)
    % plot ppc vs. spike rate
    figure(1); hold on
    plot(spikeRates, ppc(i, :), 'o-', 'Color', rand_rgb(i,:), 'MarkerFaceColor', rand_rgb(i,:))
    xlabel('spike rate (Hz)');ylabel('ppc')
    lgd1 = legend('0.0', '0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.9', '1.0');
    title(lgd1, 'spikeLockProb');
    title('Pairwise phase consistency (ppc) at 10 Hz')
    grid on

    % plot coh vs. spike rate
    figure(2); hold on
    plot(spikeRates, coh(i, :), 'o-', 'Color', rand_rgb(i,:), 'MarkerFaceColor', rand_rgb(i,:))
    xlabel('spike rate (Hz)'), ylabel('coh');
    lgd2 = legend('0.0', '0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.9', '1.0');
    title(lgd2, 'spikeLockProb');
    title('Coherence (coh) at 10 Hz')
    grid on
end

% plots 2
rand_rgb = rand(numel(spikeLockProbs),3);
for i = 1:numel(spikeLockProbs)
    % plot ppc vs. spike rate
    figure(3); hold on
    subplot(numel(spikeLockProbs), 1, i);
    plot(spikeRates, ppc(i, :), 'o-', 'MarkerEdgeColor', rand_rgb(i,:))
    xlabel('spike rate (Hz)');ylabel('ppc')
    title(sprintf('spikeLockProb = %0.2f', spikeLockProbs(i)))
    grid on

    % plot coh vs. spike rate
    figure(4); hold on
    subplot(numel(spikeLockProbs), 1, i);
    plot(spikeRates, coh(i, :), 'o-', 'MarkerEdgeColor', rand_rgb(i,:))
    xlabel('spike rate (Hz)'), ylabel('coh');
    title(sprintf('spikeLockProb = %0.2f', spikeLockProbs(i)))
    grid on
end
    