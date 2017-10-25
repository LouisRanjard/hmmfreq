function haplo = hmmfreq_dirmn2Dvq(options)
%% assemble haplotypes using HMM on windows based analysis
% 1/ Vector Quantization to align reads to reference and estimate coverage
% 2/ Multinomial
% 3/ Dirichlet multinomial from the 3 seq spots with known frequency
%
% addpath('/home/louis/Documents/Matlab/mfiles/hmmfreq');
% addpath(genpath('/home/louis/Documents/Matlab/mfiles/hmmfreq/MGLM-master'));
% addpath('/home/louis/Documents/Matlab/mfiles/nucleoveq');
%

fprintf(1,'%s - hmmfreq_dirmn2Dvq starting...\n', datestr(now)) ;

% check input
try prewinsiz=options.prewinsiz; catch prewinsiz=200; end
try ovlap=options.ovlap; catch ovlap=0.99; end
try siz=options.siz; catch siz=50; end
try rdsingle=options.rdsingle; catch rdsingle='/home/louis/Documents/Projects/Pooling3/Sequencing/Assemblies/hmm_freq/sample11/11.singleends.txt'; end
try rdpair=options.rdpair; catch rdpair='/home/louis/Documents/Projects/Pooling3/Sequencing/Assemblies/hmm_freq/sample11/11.pairends.txt'; end
try error_count=options.error_count; catch error_count=0.05; end

% load reads
reads1=fastqread(rd1);
reads2=fastqread(rd2);

% number of sub sequences (e.g. weights matrices) to assign
nsub = 3 ;

% use competitive learning to map reads to weight templates
weight_seq = alignreadsvq(reads1,reads2,'/home/louis/Documents/Projects/Pooling3/Macropodidae/EasternGrey_NC027424_Amplicons.fasta',nsub) ;
%plot(weight_seq(1).coverage);hold on;plot(weight_seq(2).coverage);plot(weight_seq(3).coverage);hold off;

% build Multiple Sequence Alignment if sequences do not share the same length
msa = multialign({weight_seq.Sequence},'terminalGapAdjust',true) ;
for n=1:nsub
  weight_seq(n).Sequence = msa(n,:) ; % copy the sequences containing gaps
  % adjust the coverage indexes with gaps
  newcoverage = zeros(1,length(msa(n,:))) ;
  p=1 ;
  for m=1:length(msa(n,:))
      if (msa(n,m)=='-')
          newcoverage(m) = nan ;
      else
          newcoverage(m) = weight_seq(n).coverage(p) ; 
          p = p+1 ;
      end
  end
  % fill the NaN positons
  weight_seq(n).coverage = fillmissing(newcoverage,'linear') ;
end

% TO DO : get the subsequences for each window with coverage
%[ cellseqs, subseqscov] = subseqset(reads,wstart,wstart+siz,nsub) ;

% memory allocation for window based analysis
lag = round(siz*(1-ovlap)) ;
if lag==0, lag=1; end % at minimum, move every position in the alignment
disp(['lag = ' num2str(lag)]) ;
xtl = 1:lag:(max(length(weight_seq(1).Sequence))-siz) ;
Y = [1 1 1;1 1 2;1 2 2;1 2 1;1 2 3] ; % combinations for the nsub most frequent ORDERED ( count(1)>count(2)>count(3) )
threehaplo = 5 ; % model for which we have exactly 3 haplotypes
y = Y ;
% find the set of rows to consider for different set size of subsequences
rowset = cell(1,nsub) ;
for n=1:nsub % get the indexes of the combination row(s) to exclude for each possible number of subsequences
    [tmprow, ~] = find(Y>n) ;
%     rowset{n} = unique(tmprow) ; % rows to exclude
    rowset{n} = setdiff(1:size(Y,1),unique(tmprow)) ; % rows to include
end
emission = NaN(size(Y,1),length(xtl)) ; % probability for the nsub^3 possible states for each window
transition = NaN(size(Y,1),length(xtl)) ; % transition between states
seqs = cell(nsub,length(xtl)) ; % save the subsequences strings
tr = zeros(1,nsub) ; % transition cost (sequence edit distance or sequence similarity pariwise) for each sequence

% expected probabilities (proportions)
true_freq = [.625 .25 .125] ; % keep the 3 most frequent subsequences
hap_count = zeros(max(length(weight_seq(1).Sequence)),3) ; % record counts of each haplotype at each position on the msa
modsel = zeros(max(length(weight_seq(1).Sequence)),2) ; % record which MODel has been SELected at each position with Likelihood

% First go through the msa to learn the local frequencies using Multinomial
% 
for wstart=xtl
    %disp([ datestr(now) ' - ' num2str(wstart) '/' num2str(xtl(end)) ])
    observations = zeros(1,nsub) ;
    [ cellseqs, subseqscov] = subseqset(weight_seq,wstart,wstart+siz,nsub) ;
    for seqsid = 1:length(subseqscov)
        seqs(seqsid,wstart) = cellseqs(seqsid) ;
        observations(seqsid) = observations(seqsid) + subseqscov(seqsid) ; % observed coverages
    end
    emissionproba = nan(1,size(Y,1)) ;
    for m=rowset{length(subseqscov)}
        model = zeros(1,length(true_freq)) ;
        unik = unique(y(m,:)) ;
        for m2 = unik
            for m1 = 1:length(true_freq)
                model(m2) = model(m2) + true_freq(m1)*(y(m,m1)==m2)  ;
            end
        end
        errora=0.01; % allow 5% error
        if (sum(model==0)==2)
            model(model==1)=1-errora;
            model(model==0)=errora/2;
        elseif (sum(model==0)==1)
            model(model>0)=model(model>0)-errora/2;
            model(model==0)=errora;
        end
        %model = sort(model,'descend') ;
        emissionproba(m) = logmnpdf(observations,model) ;
        if options.verbose==1
            fprintf(1,'%.5f %.5f %.5f\n',observations./sum(observations));
            fprintf(1,'%.5f %.5f %.5f | p = %.4f\n',model,emissionproba(m));
        end
    end
    % record the count for each haplotype
    [modsel(wstart,1),modsel(wstart,2)] = max(emissionproba) ;
    if (modsel(wstart,2)==threehaplo) % best model is 3 haplotypes, therefore we can record their respective coverages
        hap_count(wstart,1) = observations(y(modsel(wstart,2),1)) * true_freq(1) / sum(true_freq(y(modsel(wstart,2),:)==y(modsel(wstart,2),1))) ;
        hap_count(wstart,2) = observations(y(modsel(wstart,2),2)) * true_freq(2) / sum(true_freq(y(modsel(wstart,2),:)==y(modsel(wstart,2),2))) ;
        hap_count(wstart,3) = observations(y(modsel(wstart,2),3)) * true_freq(3) / sum(true_freq(y(modsel(wstart,2),:)==y(modsel(wstart,2),3))) ;
    end
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
if options.smoothing==1
    % Smoothing of coverage using FFT
    hap_countC=zeros(length(hap_count(:,1)),3);
    for m=1:3
        fft_values = fft(hap_count(:,m));
        mean_value = mean(abs(fft_values));
        threshold  = options.smoothing_factor*mean_value;
        fft_values(abs(fft_values) < threshold) = 0;
        hap_countC(:,m) = ifft(fft_values);
    end
    if options.verbose==1
        plot(hap_count(:,1));hold on;plot(hap_countC(:,1));
        plot(hap_count(:,2));plot(hap_countC(:,2)); 
        plot(hap_count(:,3));plot(hap_countC(:,3));hold off;
    end
    hap_count=hap_countC;
end

% find the maximum likelihood region to start extension from (highest confidence in coverage value of each haplotype)
modsel(modsel(:,2)~=threehaplo,1) = -Inf ; % set to -Infinity the region that do not have 3 haplotypes
[ML_multinomial,startpos] = max(modsel(:,1)) ;
% section3_s = find(diff(modsel==6)>0); % find position of set of windows where each haplotype is defined
% section3_e = find(diff(modsel==6)<0);
%pid = fopen(['./outfile' num2str(pass) '.txt'],'w') ;

% use a Gaussian Window to weight local coverage
regsiz = prewinsiz ; % size for the window to be used for measuring local coverage
%gwin = gausswin(2*siz+1) ; % initialise a Gaussian window for weighting coverage

% First, go through the msa FORWARD using hmm starting from positions with known frequencies to reconstruct haplotypes
pos_s = startpos ;
% if forward==numel(section3_e)
%     pos_e = max([reads.ends]) ;
% else
%     pos_e = pos_s + round((section3_s(forward+1)-pos_s)/2) ; % go up to middistance of gap between section3_e and the next section3_s
% end
%n=1;
for wstart=pos_s:xtl(end)
    fprintf(1,'%s - %d/%d', datestr(now), wstart, xtl(end)) ;
    % get the unique subsequences with coverage of the current window
    observations = zeros(1,nsub) ;
    [ cellseqs, subseqscov] = subseqset(weight_seq,wstart,wstart+siz,nsub) ;
    for seqsid = 1:length(subseqscov)
        seqs(seqsid,wstart) = cellseqs(seqsid) ;
        observations(seqsid) = observations(seqsid) + subseqscov(seqsid) ; % observed coverages
    end
    %observations = log(observations+1) ; % stupid
    %y = list_combi(length(subseqscov)) ; % keep the 3 most frequent or less, subseqscnt contains the number of sequences
    if wstart>1 && wstart>pos_s % number of combinations of the previous window
        %yp = list_combi(sum(cellfun(@length,seqs(:,n-1))>0)) ;
        rowsetprev = rowset{ sum(cellfun(@length,seqs(:,wstart-1))>0) } ;
    end
    if strcmp(options.prewinshape,'rectangle_prev')==1 % sum the counts for each haplotype over the REGSIZ previous windows (rectangular window)
        deb = max(1,wstart-regsiz*lag) ;
        fin = wstart ;
        region_count = sum(hap_count(deb:fin,:),1) ; % sum counts over windows
    elseif strcmp(options.prewinshape,'rectangle_around')==1 % sum the counts for each haplotype over the REGSIZ surrounding windows (rectangular window)
        deb = max(1,wstart-regsiz*lag) ;
        fin = min(wstart+regsiz*lag,xtl(end)) ;
        region_count = sum(hap_count(deb:fin,:),1) ; % sum counts over windows
    elseif strcmp(options.prewinshape,'triangle_around')==1 % apply triangular weighting window of regsiz windows before and after
        deb = max(1,wstart-regsiz) ;
        fin = min(wstart+regsiz,xtl(end)) ;
        % triangular
        region_count(1) = sum( hap_count(deb:fin,1) ./ [(wstart-deb):-1:1 1 1:(fin-wstart)]' ) ;
        region_count(2) = sum( hap_count(deb:fin,2) ./ [(wstart-deb):-1:1 1 1:(fin-wstart)]' ) ;
        region_count(3) = sum( hap_count(deb:fin,3) ./ [(wstart-deb):-1:1 1 1:(fin-wstart)]' ) ;
        % gaussian
%         region_count(1) = sum( hap_count(deb:fin,1) .* gwin((wstart-deb-regsiz+1):(fin-wstart+regsiz+1)) ) ;
%         region_count(2) = sum( hap_count(deb:fin,2) .* gwin((wstart-deb-regsiz+1):(fin-wstart+regsiz+1)) ) ;
%         region_count(3) = sum( hap_count(deb:fin,3) .* gwin((wstart-deb-regsiz+1):(fin-wstart+regsiz+1)) ) ;
    end
    region_count=round(region_count); % these are counts so better keep them integer
    for m=rowset{length(subseqscov)}
        % Dirichlet-Multinomial distribution
        counts = zeros(1,size(hap_count,2)) ;
        unik = unique(y(m,:)) ;
        for m2 = unik
            for m1 = 1:size(hap_count,2)
                counts(m2) = counts(m2) + region_count(m1)*(y(m,m1)==m2) ; % divide by sum(y(m,:)==y(m,m1))???
            end
        end
        % add errors
        %counts(counts==0) = 1 ;
        counts(counts==0) = round(error_count*sum(counts)) ; % allow 5% error in counts
        emissionproba = dirmnpdfln( observations , counts ) ;
        %fprintf(1,'Mult %.2f %.2f %.2f  => %.3f || Dirmn %.2f %.2f %.2f  => %.3f\n',model,emissionprobamn,counts,emissionproba) ;
        if wstart>1 && wstart>pos_s % percentage similarity to find transition probabilities
            %if size(yp,1)==1 % there is only one previous step, no need to compute subsequence similarities, probability is 1
            if numel(rowsetprev)==1
                I = 1 ;
                M = max(emission(:,wstart-1)) ; % all the same so all transition probabilities are equal and ==1
            else
                transitiondist = zeros(1,size(y,1))+Inf ; % initialise distance to +Infinity
                for m3=rowsetprev % consider all the possible states for the previous window
                    tr(1) = hammingd( seqs{y(m,1),wstart}, 1, siz-lag , seqs{y(m3,1),wstart-1}, 1+lag, siz ) ; % Hamming distance
                    tr(2) = hammingd( seqs{y(m,2),wstart}, 1, siz-lag , seqs{y(m3,2),wstart-1}, 1+lag, siz ) ;
                    tr(3) = hammingd( seqs{y(m,3),wstart}, 1, siz-lag , seqs{y(m3,3),wstart-1}, 1+lag, siz ) ;
                    transitiondist(m3) = sum(tr) ; % sum of mismatches = distance
                end
                %transitionproba = transitiondist'.*log(0.01) + emission(:,wstart-1) ; % using mutation probability of 0.01, proba transition = 0.01.^(distance) and take log
                transitionsim = (1+transitiondist).^(-1) ; % similarity = 1/(1+distance)
                transitionproba = log( transitionsim./sum(transitionsim(rowsetprev)) )' + emission(:,wstart-1) ;
                [M,I] = max(transitionproba) ;
            end
            transition(m,wstart) = I ;
            emission(m,wstart) = emissionproba + M ; % add emission probability
        else
            emission(m,wstart) = emissionproba ;
        end
    end
    % Save the regional coverage count: calculate the new counts for current window according to the model with highest probability
    % i.e. record the local (100 windows) coverage for each haplotype
    [ML,bestmod] = max(emission(:,wstart)) ;
    fprintf(1,' ML=%.4f\n',ML) ;
    if sum(region_count)>0
        region_freq = region_count ./ sum(region_count) ; % regional proportions (frequencies)
    else
        region_freq = [1/3 1/3 1/3] ; fprintf(1,'-- Zero Regional Coverage --\n') ;
    end
    hap_count(wstart,1) = observations(y(bestmod,1)) * region_freq(1) / sum(region_freq(y(bestmod,:)==y(bestmod,1))) ;
    hap_count(wstart,2) = observations(y(bestmod,2)) * region_freq(2) / sum(region_freq(y(bestmod,:)==y(bestmod,2))) ;
    hap_count(wstart,3) = observations(y(bestmod,3)) * region_freq(3) / sum(region_freq(y(bestmod,:)==y(bestmod,3))) ;
    %n=n+1 ;
    %fprintf(pid,'%d,%d,%d,%d,%f,%f,%f,%d\n',wstart,observations,region_count,bestmod) ;
end
%fclose(pid);

% find the best path and reconstruct the 3 haplotypes
haplo = struct('Header','','Sequence','');
haplo(1).Header = 'haplo1' ;
haplo(2).Header = 'haplo2' ;
haplo(3).Header = 'haplo3' ;
haplo(1).Sequence = repmat('N',1,max([reads.ends])) ;
haplo(2).Sequence = repmat('N',1,max([reads.ends])) ;
haplo(3).Sequence = repmat('N',1,max([reads.ends])) ;
wend = max([reads.ends]) ;
%[~,I] = nanmax(emission(:,length(xtl))) ;
[ML_forward,I] = max(emission(:,length(xtl))) ;
for x=length(xtl):-1:pos_s % index on windows
    if ~isnan(I)
        wstart = wend-lag ;
        %fprintf(1,'%d %d->%d\n',x,wstart,wend);
        if sum( cellfun(@length,seqs(:,x)) )>0 % seqs is empty when the reconstruction is partial (e.g. only forward or backward were run)
            haplo(1).Sequence(wstart:wend) = seqs{y(I,1),x}((siz-lag):siz) ;
            haplo(2).Sequence(wstart:wend) = seqs{y(I,2),x}((siz-lag):siz) ;
            haplo(3).Sequence(wstart:wend) = seqs{y(I,3),x}((siz-lag):siz) ;
        else
            fprintf(1,'%d\n',x);
        end
        wend = wstart ;
        I = transition(I,x) ;
    else
        [~,I] = max(emission(:,x-1)) ;
    end
end
if ~isnan(I) % glue the first chunk of each sequence
    haplo(1).Sequence = [ seqs{y(I,1),1}(1:(siz-1)) haplo(1).Sequence(wstart:end) ] ; % check if (siz-1) index insert a gap???
    haplo(2).Sequence = [ seqs{y(I,2),1}(1:(siz-1)) haplo(2).Sequence(wstart:end) ] ;
    haplo(3).Sequence = [ seqs{y(I,3),1}(1:(siz-1)) haplo(3).Sequence(wstart:end) ] ;
end

% Second, go through the msa BACKWARD
for wstart=pos_s:-1:xtl(1)
    fprintf(1,'%s - %d/%d', datestr(now), wstart, xtl(1)) ;
    % get the unique subsequences with coverage of the current window
    observations = zeros(1,nsub) ;
    [ cellseqs, subseqscov] = subseqset(weight_seq,wstart,wstart+siz,nsub) ;
    for seqsid = 1:length(subseqscov)
        seqs(seqsid,wstart) = cellseqs(seqsid) ;
        observations(seqsid) = observations(seqsid) + subseqscov(seqsid) ; % observed coverages
    end
    %observations = log(observations+1) ; % stupid
    %y = list_combi(length(subseqscov)) ; % keep the 3 most frequent or less, subseqscnt contains the number of sequences
    if wstart<(pos_s+1) % number of combinations of the previous window
        %yp = list_combi(sum(cellfun(@length,seqs(:,n-1))>0)) ;
        rowsetprev = rowset{ sum(cellfun(@length,seqs(:,wstart+1))>0) } ;
    end
    if strcmp(options.prewinshape,'rectangle_prev')==1 % sum the counts for each haplotype over the REGSIZ previous windows (rectangular window)
        deb = wstart ;
        fin = max(1,wstart+regsiz*lag) ;
        region_count = sum(hap_count(deb:fin,:),1) ; % sum counts over windows
    elseif strcmp(options.prewinshape,'rectangle_around')==1 % sum the counts for each haplotype over the REGSIZ surrounding windows (rectangular window)
        deb = max(wstart-regsiz*lag,1) ;
        fin = min(xtl(end),wstart+regsiz*lag) ;
        region_count = sum(hap_count(deb:fin,:),1) ; % sum counts over windows
    elseif strcmp(options.prewinshape,'triangle_around')==1 % apply triangular weighting window of regsiz windows before and after
        deb = max(1,wstart-regsiz) ;
        fin = min(wstart+regsiz,xtl(end)) ;
        % triangular
        region_count(1) = sum( hap_count(deb:fin,1) ./ [(wstart-deb):-1:1 1 1:(fin-wstart)]' ) ;
        region_count(2) = sum( hap_count(deb:fin,2) ./ [(wstart-deb):-1:1 1 1:(fin-wstart)]' ) ;
        region_count(3) = sum( hap_count(deb:fin,3) ./ [(wstart-deb):-1:1 1 1:(fin-wstart)]' ) ;
        % gaussian
%         region_count(1) = sum( hap_count(deb:fin,1) .* gwin((abs(wstart-deb-regsiz)+1):(fin-wstart+regsiz+1)) ) ;
%         region_count(2) = sum( hap_count(deb:fin,2) .* gwin((abs(wstart-deb-regsiz)+1):(fin-wstart+regsiz+1)) ) ;
%         region_count(3) = sum( hap_count(deb:fin,3) .* gwin((abs(wstart-deb-regsiz)+1):(fin-wstart+regsiz+1)) ) ;
    end
    region_count=round(region_count); % these are counts so better keep them integer
    for m=rowset{length(subseqscov)}
        % Dirichlet-Multinomial distribution
        counts = zeros(1,size(hap_count,2)) ;
        unik = unique(y(m,:)) ;
        for m2 = unik
            for m1 = 1:size(hap_count,2)
                counts(m2) = counts(m2) + region_count(m1)*(y(m,m1)==m2) ; % divide by sum(y(m,:)==y(m,m1))???
            end
        end
        % add errors
        counts(counts==0) = round(error_count*sum(counts)) ; % allow 5% error in counts
        emissionproba = dirmnpdfln( observations , counts ) ;
        %fprintf(1,'Mult %.2f %.2f %.2f  => %.3f || Dirmn %.2f %.2f %.2f  => %.3f\n',model,emissionprobamn,counts,emissionproba) ;
        if wstart<(pos_s+1) % percentage similarity to find transition probabilities
            %if size(yp,1)==1 % there is only one previous step, no need to compute subsequence similarities, probability is 1
            if numel(rowsetprev)==1
                I = 1 ;
                M = max(emission(:,wstart+1)) ; %  all the same so all transition probabilities are equal and ==1
            else
                transitiondist = zeros(1,size(y,1))+Inf ; % initialise distance to +Infinity
                for m3=rowsetprev % consider all the possible states for the previous window
                    tr(1) = hammingd( seqs{y(m,1),wstart}, 1+lag, siz, seqs{y(m3,1),wstart+1}, 1, siz-lag ) ; % Hamming distance
                    tr(2) = hammingd( seqs{y(m,2),wstart}, 1+lag, siz, seqs{y(m3,2),wstart+1}, 1, siz-lag ) ;
                    tr(3) = hammingd( seqs{y(m,3),wstart}, 1+lag, siz, seqs{y(m3,3),wstart+1}, 1, siz-lag ) ;
                    transitiondist(m3) = sum(tr) ; % sum of mismatches = distance
                end
%                 transitionproba = transitiondist'.*log(0.01) + emission(:,wstart+1) ; % using mutation probability of 0.01, proba transition = 0.01.^(distance)
                transitionsim = (1+transitiondist).^(-1) ; % similarity = 1/(1+distance)
                transitionproba = log( transitionsim./sum(transitionsim(rowsetprev)) )' + emission(:,wstart+1) ;
                [M,I] = max(transitionproba) ;
            end
            transition(m,wstart) = I ;
            emission(m,wstart) = emissionproba + M ; % add emission probability
        else
            emission(m,wstart) = emissionproba ;
        end
    end
    % Save the regional coverage count: calculate the new counts for current window according to the model with highest probability
    % i.e. record the local (100 windows) coverage for each haplotype
    [ML,bestmod] = max(emission(:,wstart)) ;
    fprintf(1,' ML=%.4f\n',ML) ;
    if sum(region_count)>0
        region_freq = region_count ./ sum(region_count) ; % regional proportions (frequencies)
    else
        region_freq = [1/3 1/3 1/3] ; fprintf(1,'-- Zero Regional Coverage --\n') ;
    end
    hap_count(wstart,1) = observations(y(bestmod,1)) * region_freq(1) / sum(region_freq(y(bestmod,:)==y(bestmod,1))) ;
    hap_count(wstart,2) = observations(y(bestmod,2)) * region_freq(2) / sum(region_freq(y(bestmod,:)==y(bestmod,2))) ;
    hap_count(wstart,3) = observations(y(bestmod,3)) * region_freq(3) / sum(region_freq(y(bestmod,:)==y(bestmod,3))) ;
    %n=n+1 ;
    %fprintf(pid,'%d,%d,%d,%d,%f,%f,%f,%d\n',wstart,observations,region_count,bestmod) ;
end
%fclose(pid);

% finish reconstructing the 3 haplotypes
[ML_backward,I] = max(emission(:,1)) ;
for x=1:pos_s % index on windows
    if ~isnan(I)
        wend = wstart+lag ;
        %fprintf(1,'%d %d->%d\n',x,wstart,wend);
        if sum( cellfun(@length,seqs(:,x)) )>0 % seqs is empty when the reconstruction is partial (e.g. only forward or backward were run)
            haplo(1).Sequence(wstart:wend) = seqs{y(I,1),x}(1:(1+lag)) ;
            haplo(2).Sequence(wstart:wend) = seqs{y(I,2),x}(1:(1+lag)) ;
            haplo(3).Sequence(wstart:wend) = seqs{y(I,3),x}(1:(1+lag)) ;
        else
            fprintf(1,'%d\n',x);
        end
        wstart = wend ;
        I = transition(I,x) ;
    else
        [~,I] = max(emission(:,x+1)) ;
    end
end
if ~isnan(I) % glue the first chunk of each sequence
    haplo(1).Sequence = [ haplo(1).Sequence(1:x-1) seqs{y(I,1),x}(1:(siz-1)) haplo(1).Sequence(x+siz:end) ] ; % check if (siz-1) index insert a gap???
    haplo(2).Sequence = [ haplo(2).Sequence(1:x-1) seqs{y(I,2),x}(1:(siz-1)) haplo(2).Sequence(x+siz:end) ] ; 
    haplo(3).Sequence = [ haplo(3).Sequence(1:x-1) seqs{y(I,3),x}(1:(siz-1)) haplo(3).Sequence(x+siz:end) ] ; 
end

% Print out the resulting sequences to fasta
fastawrite(outfile,haplo);
fprintf(1,'starting position=%d, ML_Multinomial=%.4f, ML_Forward=%.4f, ML_Backward=%.4f\n',startpos,ML_multinomial,ML_forward,ML_backward) ;

end