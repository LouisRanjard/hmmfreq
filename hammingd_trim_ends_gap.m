function [num_error]=hammingd_trim_ends_gap(seq1, seq2)
% compute Hamming distance between 2 sequences (after aligning them) while
% ignoring the gaps at the beginning and end and deleting 75bp of the first
% sequence
%
% call from command line: matlab -r "s1=fastaread('');s2=fasataread('');n=hammingd_trim_ends_gap(s1,s2);fprintf(1,'%d\n',n);quit"
%
    % align reconstructed to true, use affine gap extension penalty, trim half a read length of each reconstructed haplotype
    [~, Alignment] = nwalign(seq1.Sequence(75:end-75), seq2.Sequence, 'Alphabet', 'NT', 'GapOpen',8, 'ExtendGap',2);
    a = max( find(Alignment(1,:)~='-',1), find(Alignment(3,:)~='-',1) ) ; % ignore gaps at beginning
    b = min( find(Alignment(1,:)~='-',1,'last'), find(Alignment(3,:)~='-',1,'last') ) ; % ignore gaps at end
    % calculate Hamming distance
    num_error = hammingd( Alignment(1,:), a, b, Alignment(3,:), a, b ) ;
end