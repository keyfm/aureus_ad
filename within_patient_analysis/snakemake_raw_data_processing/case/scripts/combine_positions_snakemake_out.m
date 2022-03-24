function combine_positions_snakemake_out( path_to_list_of_input_p_files, path_to_other_p_file, path_to_output_p_file, path_to_outgroup_boolean_file, REF_GENOME_DIRECTORY, looseFQmax )
% Make a list of candiate positions

%% Version history

% % This is adapted from TDL's build_mutation_table_master_smaller_file_size_backup.m
% % Arolyn, 2018.12.19: This script was written as part of the transition
% to snakemake. It performs the part of the case step that identifies
% candidate SNP positions relative to a reference genome. 
% % Arolyn, 2019.02.21: Added input for text file indicating if sample is an
% outgroup; outgroup samples no longer input into
% generate_positions_snakemake


%% in_outgroup: booleans for whether or not each sample is in the outgroup
fprintf(1,'Processing outgroup booleans...\n');

% Input is space separated text file??
fname = [ pwd '/' path_to_outgroup_boolean_file ];
fid = fopen( fname, 'r' );
in_outgroup_string = fgetl(fid); % from echo, just one line
fclose(fid);
in_outgroup_cell = strsplit(in_outgroup_string,' ');
in_outgroup = cellfun(@(x) str2num(x), in_outgroup_cell);

fprintf(1,'Outgroup booleans:') % Print for troubleshooting
in_outgroup % Print for troubleshooting


%% Get positions on reference genome

[ChrStarts,~,~,~] = genomestats(REF_GENOME_DIRECTORY);


%% 1. Find positions with at least one fixed mutation relative to the reference genome
fprintf(1,'\n\nFinding positions with at least 1 fixed mutation...\n');

% Import list of directories for where to find variant positions for each
% sample; turn this into a cell array
fname = [ pwd '/' path_to_list_of_input_p_files ];
fid = fopen( fname, 'r' );
fstring = fgetl(fid); % from echo, just one line
paths_to_input_p_files = strsplit(fstring,' ') % currently printing out
fclose(fid);

fprintf(1,'Paths used in generate positions:') % Print for troubleshooting
paths_to_input_p_files(~in_outgroup) % Print for troubleshooting
cp = generate_positions_snakemake( paths_to_input_p_files(~in_outgroup), REF_GENOME_DIRECTORY ); % Changed 2019.02.21
%cp = generate_positions_snakemake( paths_to_input_p_files, REF_GENOME_DIRECTORY ); % Changed 2019.02.21
fprintf(1,['Found ' num2str(length(cp)) ' positions where provided vcf called a fixed variant in at least one in-group sample with FQ score < ' num2str(looseFQmax) '\n']) ;


%% 2. Find positions with within-sample polymorphisms
%fprintf(1,'\nFinding single nucleotide positions with within-sample polymorphisms...\n');
% PLACEHOLDER
fprintf(1,'\nWARNING! Finding within-sample polymorphisms has NOT been implemented in the snakemake pipeline!!!!!\n');

dp=[];
%[dp, coveragethresholds] = find_diverse_positions_no_control(loose_parameters, {SampleDirs{~in_outgroup}}, {SampleNames{~in_outgroup}}', RUN_ON_COMPUTING_CLUSTER, jobsubmitoptions_short,TEMPORARYFOLDER,SCRIPTSDIRECTORY);
%fprintf(1,'Found %i positions with within-sample polymorphism that meets loose parameters in at least 1 in-group sample \n',length(dp)) ;


%% 3. Add candidate positions manually
%fprintf(1,'\nAdding other positions to consider manually...\n');

if exist(path_to_other_p_file,'file') > 0 % check if file exists
    load(path_to_other_p_file) % variable inside mat file must be named op
    if ~exist('op','var') % error if loaded file does not contain the right variable
        error(['ERROR! File for other positions (' path_to_other_p_file ') does not contain variable op!'])
    end
else 
    op=[]; % empty vector if no positions specified manually
end
% possibly change to if statement with if nargin < N...

fprintf(1,'\nConsidering %g positions previously specified \n',length(op)) ;


%% Combine all three types of positions
fprintf(1,['\nCombining all three types of positions into one list...\n']) % Aro_Change: print this

allp = unique([dp; op; cp;]);
p = sort(allp);
p = p(p>0); %remove any 0s

%positions=p2chrpos(p,ChrStarts); % Aro_Change: not clear where this is used??


%% Save positions
fprintf(1,['\nSaving list of all positions...\n']) % Aro_Change: print this

fprintf(1,['\nSaving: ' pwd '/' path_to_output_p_file '\n']) % Aro_Change: print this
save([pwd '/' path_to_output_p_file], 'p');


