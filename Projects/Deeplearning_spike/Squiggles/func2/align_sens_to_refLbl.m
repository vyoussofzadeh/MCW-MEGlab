function sensA = align_sens_to_refLbl(sens, refLbl)
[~,ix] = ismember(refLbl, cellstr(sens.label));
sensA = sens; sensA.label = sens.label(ix); sensA.chanpos = sens.chanpos(ix,:);
end
