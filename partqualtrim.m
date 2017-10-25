function [] = partqualtrim(filefq1,filefq2,pathstr)
% addpath('/home/louis/Documents/Matlab/mfiles/hmmfreq')
% example:
% partqualtrim('/home/louis/Documents/Projects/Pooling3/Sequencing/Assemblies/nucleoveq/sample10/10_S10_L001_R1_001.fastq',...
%              '/home/louis/Documents/Projects/Pooling3/Sequencing/Assemblies/nucleoveq/sample10/10_S10_L001_R2_001.fastq',...
%              '/home/louis/Documents/Projects/Pooling3/Sequencing/Assemblies/hmm_freq/sample10/');

    % initialise seed at random
    rng('shuffle'); 

    % Load reads
    if nargin==0
        filefq1='/home/louis/Documents/Projects/Pooling3/Sequencing/Assemblies/nucleoveq/sample10/10_S10_L001_R1_001.fastq';
        filefq2='/home/louis/Documents/Projects/Pooling3/Sequencing/Assemblies/nucleoveq/sample10/10_S10_L001_R2_001.fastq';
    end
    allreads1=fastqread(filefq1);
    allreads2=fastqread(filefq2);
    
    if nargin<3
        % Get the dfirectory path, assuming both reda files are in the same directory
        pathstr = fileparts(filefq1);
    end

    % Subsample reads
    %NREAD = 13333 ; % subsample to get coverage 100x on the least abundant haplotype (assume length=5,000bp) (N*0.125*150*2)/5000=100
    for NREAD = round((5000.*[10 30 50 100 200])/(0.125*150*2))
        for n=1:10 % 10 subsets
            idr = datasample(1:numel(allreads1),NREAD,'Replace',false) ; 
            reads1=allreads1(idr);
            reads2=allreads2(idr);

            % write new fastq files
            mkdir([pathstr '/subset_' num2str(NREAD) '_' num2str(n)]) ;
            fastqwrite([pathstr '/subset_' num2str(NREAD) '_' num2str(n) '/r1.fastq'],reads1) ;
            fastqwrite([pathstr '/subset_' num2str(NREAD) '_' num2str(n) '/r2.fastq'],reads2) ;

            % Quality trimming and filtering
            cd([pathstr '/subset_' num2str(NREAD) '_' num2str(n)]);
            seqtrim([pathstr '/subset_' num2str(NREAD) '_' num2str(n) '/r1.fastq'],'Method','MaxNumberLowQualityBases','Threshold',[5 20] ) ;
            seqtrim([pathstr '/subset_' num2str(NREAD) '_' num2str(n) '/r2.fastq'],'Method','MaxNumberLowQualityBases','Threshold',[5 20] ) ;
            seqfilter([pathstr '/subset_' num2str(NREAD) '_' num2str(n) '/r1_trimmed.fastq'],'Method','MeanQuality','Threshold',20,'OutputSuffix','_filt' ) ;
            seqfilter([pathstr '/subset_' num2str(NREAD) '_' num2str(n) '/r1_trimmed_filt.fastq'],'Method','MinLength','Threshold',75,'OutputSuffix','ered' ) ;
            seqfilter([pathstr '/subset_' num2str(NREAD) '_' num2str(n) '/r2_trimmed.fastq'],'Method','MeanQuality','Threshold',20,'OutputSuffix','_filt' ) ;
            seqfilter([pathstr '/subset_' num2str(NREAD) '_' num2str(n) '/r2_trimmed_filt.fastq'],'Method','MinLength','Threshold',75,'OutputSuffix','ered' ) ;

            % use BBtools to fix the pairing between the R1 and R2 files
            system(['/home/louis/Downloads/bbmap/repair.sh '...
                'in1=' pathstr '/subset_' num2str(NREAD) '_' num2str(n) '/r1_trimmed_filtered.fastq '...
                'in2=' pathstr '/subset_' num2str(NREAD) '_' num2str(n) '/r2_trimmed_filtered.fastq '...
                'out1=' pathstr '/subset_' num2str(NREAD) '_' num2str(n) '/r1_trimmed_filtered_paired.fastq '...
                'out2=' pathstr '/subset_' num2str(NREAD) '_' num2str(n) '/r2_trimmed_filtered_paired.fastq '...
                'outsingle=' pathstr '/subset_' num2str(NREAD) '_' num2str(n) '/singletons.fastq']) ;

            % use BBtools to remove adapters
            system(['/home/louis/Downloads/bbmap/bbduk.sh '...
                'in1=' pathstr '/subset_' num2str(NREAD) '_' num2str(n) '/r1_trimmed_filtered_paired.fastq '...
                'in2=' pathstr '/subset_' num2str(NREAD) '_' num2str(n) '/r2_trimmed_filtered_paired.fastq '...
                'out1=' pathstr '/subset_' num2str(NREAD) '_' num2str(n) '/r1_adapttr_trimmed_filtered_paired.fastq '...
                'out2=' pathstr '/subset_' num2str(NREAD) '_' num2str(n) '/r2_adapttr_trimmed_filtered_paired.fastq '...
                'ref=/home/louis/Downloads/bbmap/resources/nextera.fa.gz ktrim=r k=23 mink=11 hdist=1 tpe tbo ']) ;
        end
    end
    
end
