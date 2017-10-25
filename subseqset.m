function [subseqs,cover] = subseqset(input,a,b,num,freq_limit)
% return set of subsequence (at most the 'num' most frequent) given reads and position in the alignment
% if no reads coverage in the region, return a single sequence of Ns

    if ( isfield(input,'coverage')==0 ) % reads are input
        reads = input ;
        sseqs = cell(1) ;
        wreads = reads( [reads.pos]<=a & [reads.ends]>=b ) ;
        if ~isempty(wreads)
            for p = 1:numel(wreads)
                decalade = a-wreads(p).pos ;
                sseqs{p} = wreads(p).Sequence( (1+decalade):(b-a+decalade) ) ;
            end
            [usseqs, ~, c] = unique(sseqs) ; % get unique subsequences (usseqs) and mapping from original (c)
            d = hist(c,length(usseqs)) ; % get counts
            dfreq = d./sum(d) ;          % get frequencies
            [~,idd] = sort(d,'descend') ; % sort from most frequent to least
            %num = min(num,numel(usseqs)) ; % at most, return 'num' sequences
            num = min(num,sum(dfreq>freq_limit)) ;
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
    else % alignment is input
        weight_seq = input ; 
        sseqs = cell(1,num) ;
        cov = zeros(1,num) ;
        % extract sequences and coverage of each
        for p = 1:num
            sseqs{p} = weight_seq(p).Sequence(a:b) ;
            cov(p) = sum(weight_seq(p).coverage(a:b)) ;
        end
        [usseqs, ~, c] = unique(sseqs) ; % get unique subsequences (usseqs) and mapping from original (c)
        num = min(num,numel(usseqs)) ; % update 'num' as the number of unique sequences
        d = zeros(1,num) ;
        for n=1:num % get counts
          d(n) = sum(cov(c==n)) ;
        end
        [~,idd] = sort(d,'descend') ; % sort from most frequent to least
        cover = zeros(1,num) ;
        subseqs = cell(1,num) ;
        for idx = 1:num % order 
            if (idx<=numel(d))
                cover(idx) = d(idd(idx)) ;
                subseqs{idx} = usseqs{idd(idx)} ;
            else
                break ;
            end
        end
    end
    
end