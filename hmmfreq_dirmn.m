function haplo = hmmfreq_dirmn(rdpair,rdsingle,siz,ovlap,outfile)
%% assemble haplotypes using HMM on windows based analysis
%
%addpath('/home/louis/Documents/Matlab/mfiles/hmmfreq');
%addpath(genpath('/home/louis/Documents/Matlab/mfiles/hmmfreq/MGLM-master'));


% check input
if nargin<5
    outfile=0;
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
end

% load reads
[~, n_l_s] = system(['grep -c ".$" ' rdsingle]);
nrds = str2double(n_l_s)/3 ;
if nrds>0
    rds(nrds) = struct('Header','','Sequence','','pos',0,'ends',0) ;
    fid=fopen(rdsingle);
    if fid<0
        error('Cannot open %s',rdsingle) ;
    else
    %     for rdid = 1:nrds
    %        InputText = textscan(fid,'R%*s\t%s',1) ; reads(rdid).Header =   strjoin( InputText{1} ) ;
    %        InputText = textscan(fid,'pos:\t%s',1) ; reads(rdid).pos =     str2double(InputText{1})+1 ;% indexes in input files start at 1
    %        InputText = textscan(fid,'seq:\t%s',1) ; reads(rdid).Sequence = strjoin( InputText{1} ) ;
    %     end
        InputText = textscan(fid,'R%*s\t%s\npos:\t%f\nseq:\t%s') ;
        [rds.Header] = InputText{1}{:} ;
        tmp = num2cell(InputText{2} + 1) ; [rds.pos] = tmp{:} ; % indexes in input files start at 0
        [rds.Sequence] = InputText{3}{:} ;
    end
    fclose(fid);
else
    rds = struct('Header','','Sequence','','pos',0,'ends',0) ;
end
[~, n_l_s] = system(['grep -c ".$" ' rdpair]);
nrdp = str2double(n_l_s)/3 ;
if nrdp>0
    rdp(nrdp) = struct('Header','','Sequence','','pos',0,'ends',0) ;
    fid=fopen(rdpair);
    if fid<0
        error('Cannot open %s',rdpair) ;
    else
        InputText = textscan(fid,'R%*s\t%s\npos:\t%f\nseq:\t%s') ;
        [rdp.Header] = InputText{1}{:} ;
        tmp = num2cell(InputText{2} + 1) ; [rdp.pos] = tmp{:} ; % indexes in input files start at 0
        [rdp.Sequence] = InputText{3}{:} ;
    end
    fclose(fid);
else
    rdp = struct('Header','','Sequence','','pos',0,'ends',0) ;
end
reads = [rds rdp];
% get the ends position in the alignement for each read
tmp = num2cell([reads.pos] + cellfun(@numel,{reads.Sequence}) - 1) ; [reads.ends] = tmp{:} ;


% subsample reads for tests
if 0
    NREAD = 3750 ; % subsample to get coverage 30x on the least abundant haplotype (N*0.125*150*2)/4641=30
    idr = datasample(1:numel(reads),NREAD,'Replace',false) ; 
    reads = reads(idr) ;
end

% memory allocation for window based analysis
lag = round(siz*(1-ovlap)) ;
if lag==0, lag=1; end % at minimum move every position in the alignment
disp(['lag = ' num2str(lag)]);
xtl = 1:lag:max([reads.ends]) ;
nsub = 3 ; % number of sub sequences (e.g. weights matrices) to assign
Y = list_combi(nsub) ; % get all possible combinations for the nsub most frequent
y = Y ;
% find the set of rows to consider for different set size of subsequences
rowset = cell(1,nsub) ;
for n=1:nsub % get the indexes of the combination row(s) to exclude for each possible number of subsequences
    [tmprow, ~] = find(Y>n) ;
%     rowset{n} = unique(tmprow) ; % rows to exclude
    rowset{n} = setdiff(1:size(Y,1),unique(tmprow)) ; % rows to include
end
emission = NaN(nsub^3,length(xtl)) ; % probability for the nsub^3 possible states for each window
transition = NaN(nsub^3,length(xtl)) ; % transition between states
seqs = cell(nsub,length(xtl)) ; % save the subsequences strings
tr = zeros(1,3) ; % transition cost (sequence edit distance or sequence similarity pariwise) for each sequence

% expected probabilities (proportions)
true_freq = [.625 .25 .125] ; % keep the 3 most frequent subsequences
hap_count = zeros(max([reads.ends]),3) ; % record counts of each haplotype at each position on the msa
modsel = zeros(max([reads.ends]),1) ; % record which MODel has been SELected at each position

% First go through the msa to learn the local frequencies using Multinomial
if 1
    xtl0 = 1:20:max([reads.ends]) ;
    xtl0=xtl;
    for wstart=xtl0
        disp([ datestr(now) ' - ' num2str(wstart) '/' num2str(xtl(end)) ])
        observations = zeros(1,nsub) ;
        [ cellseqs, subseqscov] = subseqset(reads,wstart,wstart+siz,nsub) ;
        for seqsid = 1:length(subseqscov)
            seqs(seqsid,n) = cellseqs(seqsid) ;
            observations(seqsid) = observations(seqsid) + subseqscov(seqsid) ; % observed coverages
        end
        emissionproba = nan(1,27) ;
        for m=rowset{length(subseqscov)}
            model = zeros(1,length(true_freq)) ;
            unik = unique(y(m,:)) ;
            for m2 = unik
                for m1 = 1:length(true_freq)
                    model(m2) = model(m2) + true_freq(m1)*(y(m,m1)==m2)  ;
                end
            end
            errora=0.01; % allow 1% error
            if (sum(model==0)==2)
                model(model==1)=1-errora;
                model(model==0)=errora/2;
            elseif (sum(model==0)==1)
                model(model>0)=model(model>0)-errora/2;
                model(model==0)=errora;
            end
            %model = sort(model,'descend') ;
            emissionproba(m) = logmnpdf(observations,model) ;
            if 0
                fprintf(1,'%.5f %.5f %.5f\n',observations./sum(observations));
                fprintf(1,'%.5f %.5f %.5f | p = %.4f\n',model,emissionproba(m));
            end
        end
        % record the count for each haplotype
        [~,modsel(wstart)] = max(emissionproba) ;
        hap_count(wstart,1) = observations(y(modsel(wstart),1)) * true_freq(1) / sum(true_freq(y(modsel(wstart),:)==y(modsel(wstart),1))) ;
        hap_count(wstart,2) = observations(y(modsel(wstart),2)) * true_freq(2) / sum(true_freq(y(modsel(wstart),:)==y(modsel(wstart),2))) ;
        hap_count(wstart,3) = observations(y(modsel(wstart),3)) * true_freq(3) / sum(true_freq(y(modsel(wstart),:)==y(modsel(wstart),3))) ;
    end
    % plot coverage total and for each individual haplotype
%     plot(sum(hap_count,2),'k'); hold on; plot(hap_count(:,1)); plot(hap_count(:,2));plot(hap_count(:,3)); hold off
    % Smoothing of coverage using root mean square
%     step=40; hap_countB=zeros(length(hap_count(:,1)),3);
%     for n=step+1:step:length(hap_count(:,1))-step
%         for m=1:3
%             hap_countB(n-step:n+step,m) = sqrt(sum(hap_count(n-step:n+step,m).^2)/sum(hap_count(n-step:n+step,m)>0)) ;
%         end
%     end
%     plot(hap_count(:,1));hold on; plot(hap_countB(:,1));
%     hap_count=hap_countB;
    % Smoothing of coverage using FFT
    hap_countC=zeros(length(hap_count(:,1)),3);
    for m=1:3
        fft_values = fft(hap_count(:,m));
        mean_value = mean(abs(fft_values));
        threshold  = 2*mean_value;
        fft_values(abs(fft_values) < threshold) = 0;
        hap_countC(:,m) = abs( ifft(fft_values) ) ; % force avoiding negative values
    end
    plot(hap_count(:,1));hold on;plot(hap_countC(:,1));
    plot(hap_count(:,2));plot(hap_countC(:,2)); 
    plot(hap_count(:,3));plot(hap_countC(:,3));hold off;
    hap_count=hap_countC;
end

for pass = 1:1 % do several pass?
    pid = fopen(['./outfile' num2str(pass) '.txt'],'w') ;
    % go through the msa to reconstruct haplotypes
    n=1;
    for wstart=xtl
        disp([ datestr(now) ' - ' num2str(wstart) '/' num2str(xtl(end)) ])
        % get the unique subsequences with coverage of the current window
        observations = zeros(1,nsub) ;
        %[ seqs(:,n), subseqscnt] = subseqset(reads,wstart,wstart+siz,nsub) ;
        [ cellseqs, subseqscov] = subseqset(reads,wstart,wstart+siz,nsub) ;
        for seqsid = 1:length(subseqscov)
            seqs(seqsid,n) = cellseqs(seqsid) ;
            observations(seqsid) = observations(seqsid) + subseqscov(seqsid) ; % observed coverages
        end
        %observations = log(observations+1) ; % stupid
        %y = list_combi(length(subseqscov)) ; % keep the 3 most frequent or less, subseqscnt contains the number of sequences
        if wstart>1  % number of combinations of the previous window
            %yp = list_combi(sum(cellfun(@length,seqs(:,n-1))>0)) ;
            rowsetprev = rowset{ sum(cellfun(@length,seqs(:,n-1))>0) } ;
        end
        % sum the counts for each haplotype over the 100 surrounding windows
        if wstart==xtl(1)
            region_count = [5 2 1] ;
        else
            deb = max(1,wstart-siz*lag) ;
            fin = min(n+siz*lag,xtl(end)) ;
            region_count = sum(hap_count(deb:fin,:),1) ; % sum counts over windows
            % apply triangular weighting window?
        end
        for m=rowset{length(subseqscov)}
    %         % Multinomial distribution on coverage for emission probabilities
    %         model = zeros(1,length(true_freq)) ;
    %         unik = unique(y(m,:)) ;
    %         for m2 = unik
    %             for m1 = 1:length(true_freq)
    %                 model(m2) = model(m2) + true_freq(m1)*(y(m,m1)==m2) ;
    %             end
    %         end
    %         %disp(model);
    %         % add errors, allow 1% errors shared
    %         errora=0.01; % allow 1% error
    %         if (sum(model==0)==2)
    %             model(model==1)=1-errora;
    %             model(model==0)=errora/2;
    %         elseif (sum(model==0)==1)
    %             model(model>0)=model(model>0)-errora/2;
    %             model(model==0)=errora;
    %         end
    %         emissionprobamn = logmnpdf(observations,model) ;
            % Dirichlet-Multinomial distribution
            counts = zeros(1,size(hap_count,2)) ;
            unik = unique(y(m,:)) ;
            for m2 = unik
                for m1 = 1:size(hap_count,2)
                    counts(m2) = counts(m2) + region_count(m1)*(y(m,m1)==m2) ; % divide by sum(y(m,:)==y(m,m1))???
                end
            end
            % add errors
            counts(counts==0) = 1 ;
            emissionproba = dirmnpdfln( observations , counts ) ;
            %fprintf(1,'Mult %.2f %.2f %.2f  => %.3f || Dirmn %.2f %.2f %.2f  => %.3f\n',model,emissionprobamn,counts,emissionproba) ;
            if wstart>1 % percentage similarity to find transition probabilities
                %if size(yp,1)==1 % there is only one previous step, no need to compute subsequence similarities, probability is 1
                if numel(rowsetprev)==1
                    I = 1 ;
                    M = 0 ; % log(1)
                else
                    %for m3=1:size(yp,1)
                    %transitiondist = zeros(1,size(yp,1)) ;
                    transitiondist = zeros(1,size(y,1))+Inf ; % initialise distance to +Infinity
                    for m3=rowsetprev % consider all the possible states for the previous window
                        tr(1) = hammingd( seqs{y(m,1),n}, 1, siz-lag , seqs{y(m3,1),n-1}, 1+lag, siz ) ; % Hamming distance
                        tr(2) = hammingd( seqs{y(m,2),n}, 1, siz-lag , seqs{y(m3,2),n-1}, 1+lag, siz ) ;
                        tr(3) = hammingd( seqs{y(m,3),n}, 1, siz-lag , seqs{y(m3,3),n-1}, 1+lag, siz ) ;
                        transitiondist(m3) = sum(tr) ; % sum of mismatches = distance
                    end
    %                 transitionsim = (1+transitiondist).^(-1) ; % similarity = 1/(1+distance)
    %                 %transitionsim = exp(-0.0001*transitiondist.^2) ; % similarity = e^(-b*D^2)
    %                 transitionproba = log( transitionsim./sum(transitionsim(rowsetprev)) )' + emission(:,n-1) ;
                    transitionproba = transitiondist'.*log(0.01) + emission(:,n-1) ; % using mutation probability of 0.01, proba transition = 0.01.^(distance)
                    [M,I] = max(transitionproba) ;
                end
                transition(m,n) = I ;
                %emission(m,n) = isnanV(emission(m,n-1)) + emissionproba + M ; % add emission probability
                emission(m,n) = emissionproba + M ; % add emission probability
            else
                emission(m,n) = emissionproba ;
            end
        end
        % Save the regional coverage count: calculate the new counts for current window according to the model with highest probability
        % i.e. record the local (100 windows) coverage for each haplotype
        [~,bestmod] = max(emission(:,n)) ;
        if sum(region_count)>0
            region_freq = region_count ./ sum(region_count) ; % regional proportions (frequencies)
        else
            region_freq = [1/3 1/3 1/3] ;
        end
        hap_count(wstart,1) = observations(y(bestmod,1)) * region_freq(1) / sum(region_freq(y(bestmod,:)==y(bestmod,1))) ;
        hap_count(wstart,2) = observations(y(bestmod,2)) * region_freq(2) / sum(region_freq(y(bestmod,:)==y(bestmod,2))) ;
        hap_count(wstart,3) = observations(y(bestmod,3)) * region_freq(3) / sum(region_freq(y(bestmod,:)==y(bestmod,3))) ;
        n=n+1 ;
        %fprintf(pid,'%d,%d,%d,%d,%f,%f,%f,%d\n',wstart,observations,region_count,bestmod) ;
    end
    fclose(pid);
end

% find the best path and reconstruct the 3 haplotypes
haplo = struct('Header','','Sequence','');
haplo(1).Header = 'haplo1' ;
haplo(2).Header = 'haplo2' ;
haplo(3).Header = 'haplo3' ;
wend = max([reads.ends]) ;
%[~,I] = nanmax(emission(:,length(xtl))) ;
[~,I] = max(emission(:,length(xtl))) ;
for x=length(xtl):-1:2 % index on windows
    wstart = wend-lag ;
    %fprintf(1,'%d %d->%d\n',x,wstart,wend);
    %y = list_combi(sum(cellfun(@length,seqs(:,x))>0)) ; % get the possible combinations for this window
    haplo(1).Sequence(wstart:wend) = seqs{y(I,1),x}((siz-lag):siz) ;
    haplo(2).Sequence(wstart:wend) = seqs{y(I,2),x}((siz-lag):siz) ;
    haplo(3).Sequence(wstart:wend) = seqs{y(I,3),x}((siz-lag):siz) ;
    wend = wstart ;
    I = transition(I,x) ;
end
% glue the first chunk of each sequence
haplo(1).Sequence = [ seqs{y(I,1),1}(1:(siz-1)) haplo(1).Sequence(wstart:end) ] ; % check if (siz-1) index insert a gap???
haplo(2).Sequence = [ seqs{y(I,2),1}(1:(siz-1)) haplo(2).Sequence(wstart:end) ] ;
haplo(3).Sequence = [ seqs{y(I,3),1}(1:(siz-1)) haplo(3).Sequence(wstart:end) ] ;

% Print out
if outfile~=0
    fastawrite(outfile,haplo);
end

end