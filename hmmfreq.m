function haplo = hmmfreq(rdpair,rdsingle,siz,ovlap,outfile)
%% assemble haplotypes using HMM on windows based analysis
% addpath('/home/louis/Documents/Matlab/mfiles/hmmfreq');


% check input
if nargin<5
    outfile='/home/louis/Documents/Projects/Pooling3/Sequencing/Assemblies/hmm_freq/sample11/S11_reconstructed_hmmf_cfullreads50.fa';
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
    %     for rdid = (nrds+1):(nrds+nrdp)
    %         InputText = textscan(fid,'R%*s\t%s',1) ; reads(rdid).Header =   strjoin( InputText{1} ) ;
    %         InputText = textscan(fid,'pos:\t%s',1) ; reads(rdid).pos =     str2double(InputText{1})+1 ;% indexes in input files start at 1
    %         InputText = textscan(fid,'seq:\t%s',1) ; reads(rdid).Sequence = strjoin( InputText{1} ) ;
    %     end
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
%seqs(:) = {repmat('N',1,siz)} ; % initialise all sequences to Ns
% temporary variables
%emissionproba = zeros(1,1) ; % emission probability for 1 window and 1 state
%transitiondist = zeros(1,nsub^3) ; % sequence hamming distance for overlap between windows
%transitionsim = zeros(1,nsub^3) ; % transform hamming distance to similarity
%transitionproba = zeros(1,nsub^3) ; % transtion probability for 1 window and 1 state from the nsub^3 previous states
tr = zeros(1,3) ; % transition cost (sequence edit distance or sequence similarity pariwise) for each sequence

% expected probabilities (proportions)
%true_freq = [.615 .24 .115 0.03] ; % keep the 4 most frequent subsequences, 3% error allowed
true_freq = [.625 .25 .125] ; % keep the 3 most frequent subsequences

% First go through the msa to learn the local frequencies
if 1
    obs_freq = zeros(xtl(end),3) ;
    obs_freq(1,:) = true_freq ; % starts at expected
    obs_freq(end,:) = true_freq ; % ends at expected
    for wstart=xtl
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
        [~,bestmod]=max(emissionproba);
        switch bestmod
            case 6  
                disp( (observations([1,2,3])./sum(observations)) ) ;
                    obs_freq(wstart,:) = (observations([1,2,3])./sum(observations)) ;
            case 8  
                disp( (observations([1,3,2])./sum(observations)) ) ;
                    obs_freq(wstart,:) = (observations([1,3,2])./sum(observations)) ;
            case 12 
                disp( (observations([2,1,3])./sum(observations)) ) ;
                    obs_freq(wstart,:) = (observations([2,1,3])./sum(observations)) ;
            case 16 
                disp( (observations([2,3,1])./sum(observations)) ) ;
                    obs_freq(wstart,:) = (observations([2,3,1])./sum(observations)) ;
            case 20 
                disp( (observations([3,1,2])./sum(observations)) ) ;
                    obs_freq(wstart,:) = (observations([3,1,2])./sum(observations)) ;
            case 22 
                disp( (observations([3,2,1])./sum(observations)) ) ;
                    obs_freq(wstart,:) = (observations([3,2,1])./sum(observations)) ;
            %otherwise
            %    if wstart>xtl(1)
            %        obs_freq(wstart,:) = obs_freq(wstart-lag,:) ;
            %    end
        end
    end
    use_freq = zeros(xtl(end),3) ;
    use_freq(1,:) = true_freq ; % starts at expected
    use_freq(end,:) = true_freq ; % ends at expected
    for n=2:(xtl(end)-1) % keep first and last values to the expected ones
        a = find(obs_freq(1:n,1)) ;
        b = find(obs_freq((n+1):end,1)) ;
        use_freq(n,1) = obs_freq(a(end),1)*(b(1)/(b(1)+n-a(end))) + obs_freq(b(1)+n,1)*((n-a(end))/(b(1)+n-a(end))) ;
        a = find(obs_freq(1:n,2)) ;
        b = find(obs_freq((n+1):end,2)) ;
        use_freq(n,2) = obs_freq(a(end),2)*(b(1)/(b(1)+n-a(end))) + obs_freq(b(1)+n,2)*((n-a(end))/(b(1)+n-a(end))) ;
        a = find(obs_freq(1:n,3)) ;
        b = find(obs_freq((n+1):end,3)) ;
        use_freq(n,3) = obs_freq(a(end),3)*(b(1)/(b(1)+n-a(end))) + obs_freq(b(1)+n,3)*((n-a(end))/(b(1)+n-a(end))) ;
    end
    if 1
      plot(use_freq(:,1)); hold on ; plot(obs_freq(:,1),'o'); ylim([0.55 0.7]); hold off;
      print([rdsingle '_amplicon1_proportion_changes.png'],'-dpng')
    end
end

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
%     for m=1:size(y,1)
    for m=rowset{length(subseqscov)}
%         if sum(m==rowset{})==0 % ignore row, 0 probability
%             emission(m,n) = -Inf ; 
%             n=n+1 ;
%             next ;
%         end
        % multinomiale on coverage for emission probabilities
        model = zeros(1,length(true_freq)) ;
        unik = unique(y(m,:)) ;
        for m2 = unik
            for m1 = 1:length(true_freq)
                %model(m2) = model(m2) + true_freq(m1)*(y(m,m1)==m2) ;
                model(m2) = model(m2) + use_freq(n,m1)*(y(m,m1)==m2) ;
            end
        end
%         model = [ true_freq(1)*(y(m,1)==1)+true_freq(2)*(y(m,2)==1)+true_freq(3)*(y(m,3)==1)...
%                   true_freq(1)*(y(m,1)==2)+true_freq(2)*(y(m,2)==2)+true_freq(3)*(y(m,3)==2)...
%                   true_freq(1)*(y(m,1)==3)+true_freq(2)*(y(m,2)==3)+true_freq(3)*(y(m,3)==3)];
        %disp(model);
        % model cannot contains errors, therefore allow 1% errors shared
        errora=0.01; % allow 1% error
        if (sum(model==0)==2)
            model(model==1)=1-errora;
            model(model==0)=errora/2;
        elseif (sum(model==0)==1)
            model(model>0)=model(model>0)-errora/2;
            model(model==0)=errora;
        end
        emissionproba = logmnpdf(observations,model) ; % Multinomial distribution
        if wstart>1e9
          fprintf(1,'%.2f %.2f %.2f  => %.3f\n',model,emissionproba) ;
        end
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
                transitionsim = (1+transitiondist).^(-1) ; % similarity = 1/(1+distance)
                %transitionsim = exp(-0.0001*transitiondist.^2) ; % similarity = e^(-b*D^2)
                %transitionproba = log( transitionsim./sum(transitionsim) )' + emission(~isnan(emission(:,n-1)),n-1) ; % turn into log probabilities and add previous path
                transitionproba = log( transitionsim./sum(transitionsim(rowsetprev)) )' + emission(:,n-1) ;
                [M,I] = max(transitionproba) ;
            end
            transition(m,n) = I ;
            %emission(m,n) = isnanV(emission(m,n-1)) + emissionproba + M ; % add emission probability
            emission(m,n) = emissionproba + M ; % add emission probability
        else
            emission(m,n) = emissionproba ;
        end
    end
    if 0
        % dynamic adaptation of true frequencies when exactly 3 subsequences are found
        uf = 0.75 ; % update factor
        [~,bestmod]=max(emission(:,n));
        switch bestmod
            case 6,  true_freq = (observations([1,2,3])./sum(observations)).*uf + [0.625 0.25 0.125].*(1-uf) ; disp(true_freq) ;
            case 8,  true_freq = (observations([1,3,2])./sum(observations)).*uf + [0.625 0.25 0.125].*(1-uf) ; disp(true_freq) ;
            case 12, true_freq = (observations([2,1,3])./sum(observations)).*uf + [0.625 0.25 0.125].*(1-uf) ; disp(true_freq) ;
            case 16, true_freq = (observations([2,3,1])./sum(observations)).*uf + [0.625 0.25 0.125].*(1-uf) ; disp(true_freq) ;
            case 20, true_freq = (observations([3,1,2])./sum(observations)).*uf + [0.625 0.25 0.125].*(1-uf) ; disp(true_freq) ;
            case 22, true_freq = (observations([3,2,1])./sum(observations)).*uf + [0.625 0.25 0.125].*(1-uf) ; disp(true_freq) ;
        end
    end
    n=n+1 ;
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
    y = list_combi(sum(cellfun(@length,seqs(:,x))>0)) ; % get the possible combinations for this window
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
% Print
fastawrite(outfile,haplo);



%% Private functions
    function [subseqs,cover] = subseqset(reads,a,b,num)
    % return set of subsequence (at most the 'num' most frequent) given reads and position in the alignment
    % if no reads coverage in the region, return a single sequence of Ns
        sseqs = cell(1) ;
        wreads = reads( [reads.pos]<=a & [reads.ends]>=b ) ;
        if ~isempty(wreads)
            for p = 1:numel(wreads)
                decalade = a-wreads(p).pos ;
                sseqs{p} = wreads(p).Sequence( (1+decalade):(b-a+decalade) ) ;
            end
            [usseqs, ~, c] = unique(sseqs) ; % get unique subsequences (usseqs) and mapping from original (c)
            num = min(num,numel(usseqs)) ; % at most, return 'num' sequences
            d = hist(c,length(usseqs)) ; % get counts
            [~,idd] = sort(d,'descend') ; % sort from most frequent to least
            cover = zeros(1,num) ;
            subseqs = cell(1,num);
            %subseqs(:) = {repmat('N',1,b-a)} ;
            for idx = 1:num
                if (idx<=numel(d))
                    cover(idx) = d(idd(idx)) ;
                    subseqs{idx} = usseqs{idd(idx)} ;
                else
                    break ;
                end
            end
        else % no coverage for this region
            subseqs = cell(1,1) ; % return a 'N' sequence
            subseqs(:) = {repmat('N',1,b-a)} ;
            cover = 0 ;
        end
    end

    function [dhamming] = hammingd(s1,a1,b1,s2,a2,b2)
        if (isempty(s1) || isempty(s2))
            dhamming = b1-a1+1 ;
        else
            dhamming = sum( (s1(a1:b1) - s2(a2:b2))~=0 ) ; % positions that are different
            N2 = (s1(a1:b1)=='N') + (s2(a2:b2)=='N') ;     % positions where each sequence is 'N'
            dhamming = dhamming + sum(N2==2) ;
        end
    end

    function [valeur] = isnanV(x)
        % function that returns the value of variable if it is not a NaN else it returns 0
        if isnan(x)
            valeur=0 ;
        else
            valeur=x;
        end
    end

end