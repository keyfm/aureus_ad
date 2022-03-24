function p = generate_positions_snakemake( paths_to_input_p_files, REF_GENOME_DIRECTORY )
% Make a list of candidate SNP positions based on variant vcf files across
% all samples


%% Version history

% % This is adapted from TDL's generate_positions.m
% % Arolyn, 2018.12.19: This script was written as part of the transition
% to snakemake. It performs the part of the case step that identifies
% candidate SNP positions relative to a reference genome. Matlab data files
% have already been created for each sample.


%% Get positions on reference genome

[ChrStarts,GenomeLength,~,~] = genomestats( REF_GENOME_DIRECTORY );


%% Get variant positions

timesvariant = zeros(GenomeLength,1); % initialize vector to count occurrances of variants across samples

%load files
for i=1:length(paths_to_input_p_files)
    %http://www.vsoch.com/2010/11/loading-dynamic-variables-in-a-static-workspace-in-matlab/
    pos=load([ pwd '/' paths_to_input_p_files{i} ]);
    if numel(pos.Positions)>2 % HC 9/13/2013
        x=chrpos2index(pos.Positions,ChrStarts);
        timesvariant(x)=timesvariant(x)+1;
    end
end

% Keep positions that vary from the reference in at least one sample but
% that don't vary from the reference in ALL samples
p = find( timesvariant > 0 & timesvariant < length(paths_to_input_p_files) );

fprintf(['Not considering ' num2str(sum(timesvariant==length(paths_to_input_p_files))) ' positions where all samples have a variant compared to the reference...\n'])

