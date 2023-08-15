function fmri_LIs = ecpfunc_read_fmri_lat()


%%
filename = '/group/jbinder/ECP/alternate/derivatives/laterality_indices/censor/all_volume_li.tsv'; % please replace with your actual path
T = readtable(filename, 'FileType', 'text', 'Delimiter', '\t');

% get unique tasks and regions
tasks = unique(T.task);
regions = unique(T.region);


li_values_all = [];
li_values_part_all = [];

% for each task and region, create a variable and assign the 'li' values
for i = 1:length(tasks)
    for j = 1:length(regions)
        task = tasks{i};
        region = regions{j};
        participant = T.participant{i};
        
        % create the variable name
        var_name = strcat(task, '_', region);
        
        idx = strcmp(T.task, task) & strcmp(T.region, region);
        % select the 'li' values for the current task and region
        li_values = T.li(idx);
        
        % assign the 'li' values to the variable
        eval([var_name ' = li_values;']);
        
        li_values_all.(var_name) = li_values;
        li_values_part_all.(var_name) = T.participant(idx);
    end
end

fmri_LIs = [];
fmri_LIs.val = li_values_all;
fmri_LIs.ID = li_values_part_all;







