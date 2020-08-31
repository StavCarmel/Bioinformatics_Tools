function CAI = calc_CAI(genes, refseq)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% This function returns the genes' CAI (codon adaptation index) of a genome.
 % genes  - the sequences of the genome's genes (ORF - open reading frame).
 % refseq - 
 %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
 %% parameters
amount_codons=64; %total 64 codons
amount_aa=21; % total 21 amino acids including the END codon

%% calculating frequency for every codon

%making a long string out of all the reference sequences
long_seq=[refseq{:}];

%frequencies
codons_freq_cell=codonbias_stav(long_seq);
codons_f=struct2cell(codons_freq_cell);

%% calculating w for each codon
w=codons_f;
for aa=1:amount_aa
    cur_f_vec=w{aa}.Freq;
    cur_max=max(cur_f_vec);
    w{aa}.Freq=cur_f_vec./cur_max;
end

%% creating a w vector and a cell of strings of the codons with the same order between them
w_vec=zeros(1,amount_codons);
cod_names=cell(1,amount_codons);
count=1;
for aa=1:amount_aa
    l=length(w{aa}.Freq);
    for c=1:l
        cur_w=w{aa}.Freq(c);
        w_vec(count)=cur_w;
        cod_names{count}=w{aa}.Codon(c);
        count=count+1;
    end
end
        
%% calculate CAI for every gene
CAI=ones(length(genes),1);

for gene = 1:length(genes)
    seq_cur=genes{gene};
    gene_length=length(genes{gene});
    w_vec_cur=zeros(1,gene_length/3);
    for nt=1:3:gene_length-2
        for codon=1:amount_codons
            if seq_cur(nt:nt+2)==cod_names{codon}{:}
               idx1=(nt+2)/3;
               w_vec_cur(idx1)= w_vec(codon);
            break
            end
        end
    end
    CAI(gene)=geomean(w_vec_cur);
end
end

