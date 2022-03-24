function generate_positions_single_sample_snakemake_withoutbool( sample_path_to_variant_vcf, sample_path_to_variant_positions, maxFQ, REF_GENOME_DIRECTORY, outgroup_boolean )

%% For debugging
fprintf(1,['\n' 'Currently examining the following vcf file: ' sample_path_to_variant_vcf '\n'])
fprintf(1,['\n' 'FQ threshold: ' num2str(maxFQ) '\n'])


%% Version History
% % Tami, 2016: An earlier version of this was not readable, re-wrote
% % Arolyn, 2018.12.19: Re-written for use with snakemake pipeline, uses
% zipped files
% % Arolyn, 2018.12.20: Deletes unzipped vcf file at the end.


%% Get reference genome information

[ChrStarts,GenomeLength,~,ScafNames] = genomestats(REF_GENOME_DIRECTORY);

% Initialize boolean vector for positions to include as candidate SNPs that
% vary from the reference genome
include = zeros(GenomeLength,1) ;


%% For outgrouop samples only:

if outgroup_boolean
    % Then save no positions
    Positions=p2chrpos(find(include),ChrStarts);
    save([pwd '/' sample_path_to_variant_positions], 'Positions');
    return % end script, no need to look at vcf file
end


%% Get candidate SNP positions from variant vcf file

% Input file
fname_in_gz=[pwd '/' sample_path_to_variant_vcf ] % Aro_Change: print zipped filename
gunzip(fname_in_gz); % Aro_Change: unzip
fname_in=fname_in_gz(1:end-3) % Aro_Change: print unzipped file name

fid=fopen(fname_in,'r');
line = fgetl(fid);

while ischar(line)
    
    if line(1)~='#'
        
        %parse line
        lt=textscan(line,'%s','Delimiter', '\t');
        l=lt{1};
        
        position_on_chr= sscanf(l{2},'%f', inf); %faster than str2double
        position=ChrStarts(strcmp(l{1},ScafNames)) + position_on_chr;
        
        alt=l{5};
        ref=l{4};
        
        %only consider for simple calls (not indel, not ambigious)
        if ~isempty(alt) & ~any(alt==',') & length(alt)==length(ref) & length(ref)==1
            
            %find and parse quality score
            xt=textscan(l{8},'%s','Delimiter', ';');
            x=xt{1};
            entrywithFQ=find(strncmp(x,'FQ',2));
            fq=sscanf(x{entrywithFQ}((find(x{entrywithFQ}=='=')+1):end),'%f', inf);  %faster than str2double
            
            if int16(fq) < maxFQ; %better than maxFQ??
                include(position)=1;
            end
        end
        
    end
    line = fgetl(fid);
    
end

fclose(fid); % close unzipped file
delete(fname_in); % delete unzipped file


%% Save positions

Positions=p2chrpos(find(include),ChrStarts);
save([pwd '/' sample_path_to_variant_positions], 'Positions');

end

