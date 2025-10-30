% classify brain regions based on MONTAGE assignments
function area = classifyArea_MAZE(areaStr)

area.isFrontal = false;
area.isMTL = false;
area.isTG = false;
area.isP = false;
area.isOc = false;
area.isSP = false;
area.isLeft = false;
area.isAm = false;

%Mtl sub-areas
area.isHip = false;
area.isEC = false;
area.isPHG = false;

if strcmpi(areaStr(1),'L')
    area.isLeft = true;
end

MTL_sub_areas = {'REC','LEC','RAH','LAH','LMH','RMH','LPHG','RPHG','LOF-AC'};
Hip_sub_areas = {'RAH','LAH','LMH','RMH','LH','RH','RHH','LHH','RTH','LTH'};
EC_sub_areas = {'REC','LEC'};
PHG_sub_areas = {'LPHG','RPHG','LOF-AC'}; 
                               

frontal_areas = {'LOF','ROF','LAF','RAF','RIF-dAC','RFA','LOPR','LAC','LOFC','ACIN','LMMFG',...
    'LPMFG','LAMFG','LMFG1','LMFG2','RMFG1','RMFG2','ROFC','RPORB',...
    'AntCingulate','AntInsula','LpostInsu','SupFrontal'};

% Can we merge these to one of the above groups?
insula = {'RLI','LLI','LTPPI','LOPPI','AINS','PINS'};

% frontal_sub_areas = {'Pz','Ez','C3','C4'};
temporal_gyrus = {'LSTG','RSTG','RaSTG','LaSTG','RpSTG','RpMTG','LpSTG','LpMTG','RTP', 'LSTG',...
                    'LHSG','LPRC', 'LSTGP'};
                                
amigdala = {'RA','LA','RBAA','LAM','RAM'};

parietal = {'RPP','RPS','RPT','RAP','LAP','LTP'};

occipital = {'RO','LO','RTO','RIO','RPO','RMO','RSO','RIOP','RTOP'};

spindleAreas = {'RPP','RPS','RPT','RAP','LAP','LTP','RAC','LAC','LPC','RPC','RpSMA','LpSMA',...
                'RO','LO','RTO','RIO','RPO','RMO','RSO','RIOP',...
                'RPC','Pz','C4'};


isMTL = 0;isFrontal = 0; isTG = 0; isSP = 0;
for ii = 1:length(MTL_sub_areas)
    if strfind(areaStr,MTL_sub_areas{ii})
        area.isMTL = true;
    end
end


% Frontal areas - include insula
for ii = 1:length(frontal_areas)
    if strcmpi(areaStr,frontal_areas{ii}) 
        area.isFrontal = true;
    end
end
for ii = 1:length(insula)
    if strcmpi(areaStr,insula{ii}) 
        area.isFrontal = true;
    end
end

% Temporal areas - temporal gyrus, amigdala
for ii = 1:length(temporal_gyrus)
    if strcmpi(areaStr,temporal_gyrus{ii})
        area.isTG = true;
    end
    
end
% for ii = 1:length(amigdala)
%     if strcmpi(areaStr,amigdala{ii})
%         area.isTG = true;
%     end
%     
% end

for ii = 1:length(amigdala)
    if strcmpi(areaStr,amigdala{ii})
        area.isAm = true;
        area.isMTL = true;
    end    
end

for ii = 1:length(spindleAreas)
    if strcmpi(areaStr,spindleAreas{ii})
        area.isSP = true;
    end
end

for ii = 1:length(parietal)
    if strcmpi(areaStr,parietal{ii})
        area.isP = true;
    end
end

for ii = 1:length(occipital)
    if strcmpi(areaStr,occipital{ii})
        area.isOc = true;
    end
end

for ii = 1:length(Hip_sub_areas)
    if strcmpi(areaStr,Hip_sub_areas{ii})
        area.isHip = true;
        area.isMTL = true;
    end
end

for ii = 1:length(EC_sub_areas)
    if strcmpi(areaStr,EC_sub_areas{ii})
        area.isEC = true;
        area.isMTL = true;
    end
end

for ii = 1:length(PHG_sub_areas)
    if strcmpi(areaStr,PHG_sub_areas{ii})
        area.isPHG = true;
    end
end


end