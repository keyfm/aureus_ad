function pileup_to_diversity_matrix_snakemake( sample_path_to_pileup, sample_path_to_diversity, sample_path_to_coverage, REF_GENOME_DIRECTORY )

%% THIS VERSION IS INTENDED FOR CLONES AND NOT METAGENOMIC SNP FINDING

%% Version History
% Arolyn created on 2018.12.18 from pileup_to_diversity_matrix
% Modified to work with snakemake!!
% Tami modified 2019.01 to change behavior when indels are found -->
% ...record call information at location... this new behavior is NOT good
% ...for metagenomic analysis (intrasample polymorphisms)
% Tami 2019.01 REMOVED calculation of statistics for metagenomics/intrasample
% SNPs, but left blanks for compatability with previous versions
% Tami 2019.01 ADDED calculation of just deletions, makes it backward
% compatible while also being able to track deletions and insertions
% seperately

%% Some notes

%fname_in should me the name of a pileup file output from mpileup, including
%base qualities, mapping qualities, and tail distances reported. The
%parents of this script (pileup_to_diversity_matrix) used this information

%fname_out is what you want to save the resulting matrix as, I call it 'data.mat'

%ChrStarts is an array holding the indices in the position dimension
%corresponding to the start of a new chromsome. For B. dolosa this is [0,3413074,5589203]

%[A T C G a t c g Aq ... gq Am .... gm  At .... gt Ps Pb Pm Pftd Prtd E D]
% Aq is the average phred qualities of all A's
% Am is the average mapping qualities of all A's
% At is the average tail distance of all A's
% Ps is the p value for strand bias (fishers test)
% Pb is the p value for the base qualities being the same for the two
% different types of calls (1st major, 2nd major nt, either strand) (ttest)
% Pm is the p value for the mapping qualities being the same for the two
% different types of calls (ttest)
% Pftd is the p value for the tail distantces on the forward strand
% being the same for the two different types of calls (ttest)
% Pftd is the p value for the tail distantces on the reverse strand
%being the same for the two  different types of calls (ttest)
% E is number of calls at ends of a read
% D is number of reads supporting indels in the +/- (indelregion) bp region


%% Important constants
Phred_offset=33; %mpileup is told that the reads are in fastq format and corrects that so its always at Phred+33 when the mpileup comes out
nts=int16('ATCGatcg');
indels=int16('+-');
numfields=40; %CHANGED 2019.01
numtests=50; %how many ttests to do in parallel
indelregion=3;


%% Input and output file names

% Path to strain.pileup (input)
fname_in=[pwd '/' sample_path_to_pileup ] % Aro_Change: new; print this
% Path to diversity.mat (output)
fname_out=[pwd '/' sample_path_to_diversity ] % Aro_Change: new; print this
% Path to coverage.mat (secondary output)
fname_out_cov=[pwd '/' sample_path_to_coverage ]; % Aro_Change: new; print this
%fname_in='strain.pileup'; % Aro_Change: old
%fname_out='diversity.mat'; % Aro_Change: old


%% Get reference genome information

%load alignment_info % Aro_Change: old
%refgenome=ae.Genome; % Aro_Change: old

fprintf(1,['\n Reference genome: ' REF_GENOME_DIRECTORY '\n'])
[ChrStarts,GenomeLength,ChromosomeIndicator,ScafNames] = genomestats(REF_GENOME_DIRECTORY); % Aro_Change: new
fprintf(1,['\n' num2str(GenomeLength) '\n']) % Aro_Change: print this
%[ChrStarts,GenomeLength,ChromosomeIndicator,ScafNames]=genomestats(['/scratch/mit_lieberman/Reference_Genomes/'
%refgenome]); %RUN_ON_CLUSTER % Aro_Change: old


%% Initialize

data=zeros(numfields,GenomeLength,'uint16'); %[A T C G a t c g Aq ... gq Am .... gm  At .... gt Ps Pb Pm Pftd Prtd E D]
temp=zeros(numfields,1); %holds data for each line before storing in data

fid=fopen(fname_in);

line = fgetl(fid);
linec = 0; 

while ischar(line)
	%fprintf(1, [line '\n']);     
    %parse line
    lt=textscan(line,'%s','Delimiter', '\t'); l=lt{1};
    chr=l{1}; 
    
    if numel(ChrStarts)==1
        position= sscanf(l{2},'%f', inf);
    else
        if sum(ismember(ScafNames,chr))==0
	    ScafNames'
	    chr
            error('Scaffold name in pileup file not found in reference')
        end
        %e.g. define position by looking at 28th character in string defining B. dolosa chrosome
        position= double(ChrStarts(ismember(ScafNames,chr))+sscanf(l{2},'%f'));
    end
  
    
    ref=find(nts==int16(l{3}));
    if ref > 4
        ref = ref - 4;
    end

    calls=int16(l{5});
    bq=int16(l{6}); %base quality, BAQ corrected
    mq=int16(l{7}); %mapping quality
    td=sscanf(l{8},'%f,'); %distance from tail
    
    %starts of reads
    startsk=find(calls==94); %'^'
    for k=startsk
        calls(k:k+1)=-1; %remove mapping character, absolutely required because the next chracter could be $
    end
    
    %read ends
    endsk=find(calls==36); %'$'
    temp(end-1)=length(endsk)+length(startsk); %record how many calls were at the start/end of a read
    calls(endsk)=-1;
    
    %indels and calls from reads supporting indels
    indelk=[find(calls==43) find(calls==45)]; %'-+'
    for k=indelk
        %find the size of the indel
        if calls(k+2) >= 48 && calls(k+2) < 58 %indel of size >9 and <100
            indelsize=str2double(char(calls(k+1:k+2)));
            indeld=2;
        else %indel size <9
            indelsize=double(calls(k+1)-48);
            indeld=1;
        end
        %record that indel was found in nearby region
        %store directly into data, as it affects lines earlier and later
        if calls(k)==45 %-deletion
            if (position-indelregion+1) > 0 && (position+indelsize+indelregion) <= GenomeLength
                data(39:40,position-indelregion+1:position+indelsize+indelregion)=data(39:40,position-indelregion+1:position+indelsize+indelregion)+1;
            elseif (position-indelregion+1) > 0
                data(39:40,position-indelregion+1:end)=data(39:40,position-indelregion+1:end)+1;
            else
                data(39:40,1:position+indelsize+indelregion)=data(39:40,1:position+indelsize+indelregion)+1;
            end
        else %insertion
            %insertion isn't indexed on chromosome, no need for anything
            %complicated
            if (position-indelregion+1) > 0 && (position+indelregion) <= GenomeLength
                data(39,position-indelregion+1:position+indelregion)=data(39,position-indelregion+1:position+indelregion)+1;
            elseif (position-indelregion+1) > 0
                data(39,position-indelregion+1:end)=data(39,position-indelregion+1:end)+1;
            else
                data(39,1:position+indelsize)=data(39,1:position+indelsize)+1;
            end
        end
        %remove calls from counting
        calls(k:(k+indeld+indelsize))=-1; %don't remove base that precedes an indel
       % COMMENTED OUT 2019.01 TDL calls(k-1)=300; %keep this number for indexing properly, but don't count it towards ATCG
    end
    %ignore '*', as these deletion markers won't be counted towards score
    %qualities, etc and the presence of a deletion is recorded by the
    %upstream '-' indicator
    
    
    
    %replace all reference with their actual calls
    if ref
        calls(calls==46)=nts(ref); %'.'
        %(1,[num2str(ref) ',' num2str(position) '\n']); %for error tracking
        calls(calls==44)=nts(ref+4); %','
    end
    
    %index reads for finding scores
    simplecalls=calls((calls>0)); %simplecalls is a transformation of calls such that each calls position in simplecalls corresponds to its position in bq, mq, td
    
    
    
    %calculate how many of each type and average scores
    for nt=1:8
        if ~isempty(find(simplecalls==nts(nt),1))
        %    fprintf(1,[num2str(length(simplecalls)) ',' num2str(temp(nt)) ',' num2str(length(nts)) ',' num2str(length(bq)) ',' num2str(length(mq)) ',' num2str(length(td))]);
                        
            temp(nt)=sum(simplecalls==nts(nt));
            temp(nt+8)= int16(sum(bq(simplecalls==nts(nt)))/temp(nt))- Phred_offset;
            temp(nt+16)=int16(sum(mq(simplecalls==nts(nt)))/temp(nt)) - 33;
            temp(nt+24)=int16(sum(td(simplecalls==nts(nt)))/temp(nt));
        end
    end
    

%% This block commented out because it is not needed for simple single-colony samples 
%    

    %find major and nextmajor allele
    [~, sortedpositions] = sort(temp(1:4)+temp(5:8));
    n1 = sortedpositions(end);
    n2 = sortedpositions(end-1);
    
%calculate stats which require distributions during loop, strand bias
%     %can de done late
%     x=(simplecalls==nts(n1) | simplecalls==nts(n1+4));
%     y=(simplecalls==nts(n2) | simplecalls==nts(n2+4));
%     if (sum(temp(1:4)) > 20 && sum(temp(5:8))> 20 && sum(y)*200>sum(x))  %only calcualte p values of there are greater than 20 reads on each strand and MAF < .995
%         
%         bttests_x=bq(x)-Phred_offset; bttests_y=bq(y)-Phred_offset;
%         mttests_x=mq(x)-Phred_offset; mttests_y=mq(y)-Phred_offset;
%         fttests_x=td(simplecalls==nts(n1)); fttests_y=td(simplecalls==nts(n2));
%         rttests_x=td(simplecalls==nts(n1+4)); rttests_y=td(simplecalls==nts(n2+4));
%         
%         [~,bp]=ttest2(double(bttests_x),double(bttests_y));
%         [~,mp]=ttest2(double(mttests_x),double(mttests_y));
%         [~,fp]=ttest2(double(fttests_x),double(fttests_y));
%         [~,rp]=ttest2(double(rttests_x),double(rttests_y));
%         
%         temp(end-5)=int16(-log10(bp));
%         temp(end-4)=int16(-log10(mp));
%         temp(end-3)=int16(-log10(fp));
%         temp(end-2)=int16(-log10(rp));
%         
%         %strand bias
%         p=fexact([temp(n1) temp(n2); temp(n1+4) temp(n2+4)], 0); %strand bias
%         temp(end-6)=int16(-log10(p(3)));
%         
%         
%         
%     end
 %%   
    %store the data!
    data(1:(end-1),position)=int16(temp(1:end-1)); %last place is for indels, can be modified by earlier and later positions
    
    
    line = fgetl(fid);
    temp=zeros(numfields,1);
	linec = linec+1;     
    
end

fclose(fid);


coverage=sum(data(1:8,:));

% Save coverage.mat
save(fname_out_cov,'coverage');

% Save diversity.mat
save(fname_out,'data');

return

