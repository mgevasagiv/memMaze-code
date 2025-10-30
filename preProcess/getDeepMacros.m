% getDeepMacros - a list of medial contacts for a specific pt based on
% macro-montage.
% input: 
% removeBad removes both high spike activity and bad channels.
%
% %% Channel flag codes
% % badChannel:
%     1. General badness (can be either crazy activity or too many NaNs, p489 had 6),
%     2. Large reactions to stimulations (p489 didn't have such channels, unless I missed them),
%     3. Nearly identical to next contact (correlation > 0.99),
%     4. Opposite polarity to the other contacts in the electrode.
% highSpikesActivity = high spikes activity which I marked in a different field
% 
% Related: see getChannelListMacros() for all contacts
function channelList = getDeepMacros(MacroMontage, AreaStr,removeBad)


if nargin < 3
    removeBad = 1;
end

if nargin < 2
    getAllChan = 1;
else
    getAllChan = 0;
end

if getAllChan
    channelList= [];
    for ii = 1:length(MacroMontage)
        if (    strcmpi(MacroMontage(ii).Area,'Pz') || ...
                strcmpi(MacroMontage(ii).Area,'Cz') || ...
                strcmpi(MacroMontage(ii).Area,'Ez') || ...
                strcmpi(MacroMontage(ii).Area,'C3') || ...
                strcmpi(MacroMontage(ii).Area,'C4')  ); break; end;
        
        if ii == 1 | ~strcmp(MacroMontage(ii).Area,MacroMontage(ii-1).Area)
            if ( ~strcmpi(MacroMontage(ii).Area,'NaN') )
                if(isfield(MacroMontage(ii),'badChannel') && ~isempty(MacroMontage(ii).badChannel))
                    if isempty(MacroMontage(ii+1).badChannel)
                        channelList = [channelList, ii+1];
                    elseif ~isempty(MacroMontage(ii+2).badChannel)
                        channelList = [channelList, ii+2];
                    elseif ~isempty(MacroMontage(ii+3).badChannel)
                        channelList = [channelList, ii+3];
                    else
                        warning('4 consecutive bad channels - skipping bad electrode')
                        continue
                    end
                else
                    channelList = [channelList, ii];
                end
                
            end
        end
    end
    
else
    
    channelList= [];
    for ii = 1:length(MacroMontage)
        if (    strcmpi(MacroMontage(ii).Area,'Pz') || ...
                strcmpi(MacroMontage(ii).Area,'Cz') || ...
                strcmpi(MacroMontage(ii).Area,'Ez') || ...
                strcmpi(MacroMontage(ii).Area,'C3') || ...
                strcmpi(MacroMontage(ii).Area,'C4')  ); break; end;
        
        if ii == 1 | ~strcmp(MacroMontage(ii).Area,MacroMontage(ii-1).Area)
            
            
            if ( ~strcmpi(MacroMontage(ii).Area,'NaN') )
                
                if strcmpi(AreaStr,MacroMontage(ii).Area)
                    
                    if removeBad
                        
                        if(isfield(MacroMontage(ii),'badChannel') && ~isempty(MacroMontage(ii).badChannel))
                            if isempty(MacroMontage(ii+1).badChannel)
                                channelList = [ii+1];
                            elseif isempty(MacroMontage(ii+2).badChannel)
                                channelList = [ii+2];
                            else
                                warning('three consecutive bad channels - skipping bad electrode')
                                continue
                            end
                        else
                            channelList = [channelList, ii];
                        end
                        
                    else
                        channelList = [channelList, ii];
                        
                    end
                    
                end
            end
        end
    end
end
