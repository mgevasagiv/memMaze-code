load(runData(iPatient).macroMontageFileName,'MacroMontage')

for ii = 1:length(MacroMontage)
    
    channel(ii) = MacroMontage(ii).Channel;
    Area{ii} = MacroMontage(ii).Area;
    
    highSpikesActivity(ii) = ~isempty(MacroMontage(ii).highSpikesActivity);
    bad_chan(ii) = ~isempty(MacroMontage(ii).bad_chan);
end
varnames = {'channel','area','highSpikesActivity','bad_chan'};
TABLE = table(channel',Area', highSpikesActivity', bad_chan','VariableNames',varnames);

[folder filename] = fileparts(runData(iPatient).macroMontageFileName)
writetable(TABLE, fullfile(folder,sprintf('%s.csv',filename)))

