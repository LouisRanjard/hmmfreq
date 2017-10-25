function haplo = hmmfreq(rdpair,rdsingle,siz,ovlap)
%% assemble haplotypes using HMM on windows based analysis
% addpath('/home/louis/Documents/Matlab/mfiles/hmmfreq');


% check input
if nargin<4
    ovlap=0.9;
    if nargin<3
        siz=100;
        if nargin<2
            rdsingle='/home/louis/Documents/Projects/Pooling3/Sequencing/Assemblies/hmm_freq/sample11/11.singleends.txt' ;
            if nargin<1
                rdpair='/home/louis/Documents/Projects/Pooling3/Sequencing/Assemblies/hmm_freq/sample11/11.pairends.txt' ;
            end
        end
    end
end

% load reads
[~, n_l_s] = system(['grep -c ".$" ' rdsingle]); nrds = str2double(n_l_s)/3 ;
[~, n_l_s] = system(['grep -c ".$" ' rdpair]);   nrdp = str2double(n_l_s)/3 ;
N = nrds + nrdp ;
reads(N) = struct('Header','','Sequence','','pos',0,'ends',0) ;
fid=fopen(rdsingle);
if fid<0
    error('Cannot open %s',rdsingle) ;
else
    for rdid = 1:nrds
        InputText = textscan(fid,'R%*s\t%s',1) ; reads(rdid).Header =   strjoin( InputText{1} ) ;
        InputText = textscan(fid,'pos:\t%s',1) ; reads(rdid).pos =     str2double(InputText{1})+1 ;% indexes in input files start at 1
        InputText = textscan(fid,'seq:\t%s',1) ; reads(rdid).Sequence = strjoin( InputText{1} ) ;
    end
end
fclose(fid);
fid=fopen(rdpair);
if fid<0
    error('Cannot open %s',rdpair) ;
else
    for rdid = (nrds+1):(nrds+nrdp)
        InputText = textscan(fid,'R%*s\t%s',1) ; reads(rdid).Header =   strjoin( InputText{1} ) ;
        InputText = textscan(fid,'pos:\t%s',1) ; reads(rdid).pos =     str2double(InputText{1})+1 ;% indexes in input files start at 1
        InputText = textscan(fid,'seq:\t%s',1) ; reads(rdid).Sequence = strjoin( InputText{1} ) ;
    end
end
fclose(fid);

 
% subsample reads for tests
if 0
    NREAD = 3750 ; % subsample to get coverage 30x on the least abundant haplotype (N*0.125*150*2)/4641=30
    idr = datasample(1:numel(reads),NREAD,'Replace',false) ; 
    reads = reads(idr) ;
end

% get the ends poistion in the alignement for each read
for r=1:numel(reads)
    reads(r).ends = reads(r).pos + length(reads(r).Sequence) - 1 ;
end

% memory allocation for window based analysis
lag = round(siz*(1-ovlap)) ;
xtl = 1:lag:max([reads.ends]) ;
nsub = 3 ; % number of sub sequences (e.g. weights matrices) to assign
y = list_combi(nsub) ; % keep the 3 most frequent
emission = zeros(nsub^3,length(xtl)) ; % probability for the nsub^3 possible states for each window
transition = zeros(nsub^3,length(xtl)) ; % transition between states
seqs = cell(3,length(xtl)) ; % save the subsequences strings
seqs(:) = {repmat('N',1,siz)} ;% initialise all sequences to Ns
% temporary variables
%emissionproba = zeros(1,1) ; % emission probability for 1 window and 1 state
transitionproba = zeros(1,nsub^3) ; % transtion probability for 1 window and 1 state from the nsub^3 previous states
tr = zeros(1,3) ; % transition cost (sequence edit distance or sequence similarity pariwise) for each sequence

% calculate probabilities
%true_freq = [.615 .24 .115 0.03] ; % keep the 4 most frequent subsequences, 3% error allowed
true_freq = [.625 .25 .125] ; % keep the 3 most frequent subsequences
n=1;
for wstart=xtl
    % get the unique subsequences with coverage of the current window
    [ seqs(:,n), subseqscnt] = subseqset(reads,wstart,wstart+siz,3) ;
    observations = subseqscnt ; % observed coverages
    for m=1:size(y,1)
        % multinomiale on coverage for emission probabilities
        % = [ sum(covers(y(m,:)==1)) sum(covers(y(m,:)==2)) sum(covers(y(m,:)==3)) ] ;
        model = [ true_freq(1)*(y(m,1)==1)+true_freq(2)*(y(m,2)==1)+true_freq(3)*(y(m,3)==1)...
                  true_freq(1)*(y(m,1)==2)+true_freq(2)*(y(m,2)==2)+true_freq(3)*(y(m,3)==2)...
                  true_freq(1)*(y(m,1)==3)+true_freq(2)*(y(m,2)==3)+true_freq(3)*(y(m,3)==3)];
        % model cannot contains errors, therefore allow 1% errors shared
        errora=0.05;
        if (sum(model==0)==2)
            model(model==1)=1-errora;
            model(model==0)=errora/2;
        elseif (sum(model==0)==1)
            model(model>0)=model(model>0)-errora/2;
            model(model==0)=errora;
        end
        emissionproba = log( mnpdf(observations,model) ) ;
        %fprintf(1,'%.2f %.2f %.2f  => %.3f\n',model,emission(m,n)) ;
        if wstart>1 % percentage similarity to find transition probabilities
            for m2=1:size(y,1) % consider all the possible states for the previous window
                tr(1) = hammingd( seqs{y(m,1),n}, 1, siz-lag , seqs{y(m2,1),n-1}, 1+lag, siz ) ; % hamming distance
                tr(2) = hammingd( seqs{y(m,2),n}, 1, siz-lag , seqs{y(m2,2),n-1}, 1+lag, siz ) ;
                tr(3) = hammingd( seqs{y(m,3),n}, 1, siz-lag , seqs{y(m2,3),n-1}, 1+lag, siz ) ;
                transitionproba(m2) = sum(tr) ;
            end
            transitionproba = log( transitionproba./sum(transitionproba) ) ; % turn into probabilities
            [M,I] = max(transitionproba) ;
            transition(m,n) = I ;
            emission(m,n) = emission(m,n-1) + M ;
        else
            emission(m,n) = emissionproba ;
        end
    end
    n=n+1;
end
% find the best path and reconstruct the 3 haplotypes
haplo = struct('Header','','Sequence','');
haplo(1).Header = 'haplo1' ;
haplo(2).Header = 'haplo2' ;
haplo(3).Header = 'haplo3' ;
wend = max([reads.ends]) ;
for x=length(xtl):-1:2 % index on windows
    %fprintf(1,'%d\n',x);
    [~,I] = max(emission(:,x)) ;
    wstart = wend-lag ;
    haplo(1).Sequence(wstart:wend) = seqs{y(I,1),x}((siz-lag):siz) ;
    haplo(2).Sequence(wstart:wend) = seqs{y(I,2),x}((siz-lag):siz) ;
    haplo(3).Sequence(wstart:wend) = seqs{y(I,3),x}((siz-lag):siz) ;
    wend = wstart ;
end
% glue the first chunk of each sequence
haplo(1).Sequence = [ seqs{y(I,1),1}(1:(siz-1)) haplo(1).Sequence(2:end) ] ;
haplo(2).Sequence = [ seqs{y(I,2),1}(1:(siz-1)) haplo(2).Sequence(2:end) ] ;
haplo(3).Sequence = [ seqs{y(I,3),1}(1:(siz-1)) haplo(3).Sequence(2:end) ] ;
% Print
fastawrite('/home/louis/Documents/Projects/Pooling3/Sequencing/Assemblies/hmm_freq/sample11/S11_reconstructed_hmmf_b.fa',haplo);



%% Private functions
    function [subseqs,cover] = subseqset(reads,a,b,num)
    % return set of subsequence (N most frequent) given reads and position in the alignment
        sseqs = cell(1) ;
        wreads = reads( [reads.pos]<=a & [reads.ends]>=b ) ;
        if ~isempty(wreads)
            for p = 1:numel(wreads)
                decalade = a-wreads(p).pos ;
                sseqs{p} = wreads(p).Sequence( (1+decalade):(b-a+decalade) ) ;
            end
            [usseqs, ~, c] = unique(sseqs) ; % get unique subsequences (usseqs) and mapping from original (c)
            d = hist(c,length(usseqs)) ; % get counts
            [~,idd] = sort(d,'descend') ; % sort from most frequent to least
            cover = zeros(1,num) ;
            subseqs = cell(1,num);
            subseqs(:) = {repmat('N',1,b-a)} ;
            for idx = 1:num
                if (idx<=numel(d))
                    cover(idx) = d(idd(idx)) ;
                    subseqs{idx} = usseqs{c(idd(idx))} ;
                else
                    break ;
                end
            end
        else % no coverage for this region
            subseqs = cell(1,num) ;
            subseqs(:) = {repmat('N',1,b-a)} ;
            cover = zeros(1,num) ;
        end
    end

    function [dhamming] = hammingd(s1,a1,b1,s2,a2,b2)
        if (isempty(s1) || isempty(s2))
            dhamming = b1-a1 ;
        else
            dhamming = sum( s1(a1:b1) - s2(a2:b2) ) ;
        end
    end


end