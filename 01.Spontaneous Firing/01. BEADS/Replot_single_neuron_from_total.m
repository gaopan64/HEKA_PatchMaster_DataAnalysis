neuron_name = '';

% Intrinsic properties
sample_rate_kHz = 40;
T = 1/sample_rate_kHz;

for theFiles_i = 1 : length(theFiles)
    if strcmp(theFiles(theFiles_i).name,neuron_name) 
        neuron_index = theFiles_i;
        fprintf(1, "Found the neuron '%s' in the collection\n", neuron_name);
    else
        fprintf(1, "Search next neuron...\n");
    end
end

if exist('neuron_index','var') == 0
    fprintf(1, "Can't find the neuron '%s' in 'theFiles'.\n", neuron_name);
    return
end

tStart = tic;

baseline_cor_v = FiringCollect(:, neuron_index);

fprintf(1, "Start analyzing...\n");

SponFiring_TotalAnal_fun(baseline_cor_v,T,neuron_name)

fprintf(1, 'The whole process takes %ss\n', num2str(toc(tStart)));
