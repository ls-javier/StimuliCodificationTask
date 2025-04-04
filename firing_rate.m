function [fRates] = firing_rate(trials_spikes,timeSamples,window_size)
%FIRING_RATE Summary of this function goes here
%   Detailed explanation goes here

nTimes=length(timeSamples);
[ntrials,TTrial]=size(trials_spikes);
fRates=zeros(ntrials,nTimes);
for i=1:ntrials
    for j=1:nTimes
        ini=timeSamples(j);
        index=trials_spikes(i,ini-window_size:ini);
        index_spikes=find(index==1);
        num_spikes=length(index_spikes);
        fRates(i,j)=(num_spikes)/(window_size);
    end
end
end

