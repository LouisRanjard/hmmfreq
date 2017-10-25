function [num_error]=d_to_true(haplo, trueh)
  num_error=zeros(1,length(haplo));
  for n=1:length(haplo)
%     % align reconstructed to true, use affine gap extension penalty, trim half a read length of each reconstructed haplotype
%     [~, Alignment] = nwalign(haplo(n).Sequence(75:end-75), trueh(n).Sequence, 'Alphabet', 'NT', 'GapOpen',8, 'ExtendGap',2);
%     a = max( find(Alignment(1,:)~='-',1), find(Alignment(3,:)~='-',1) ) ; % ignore gaps at beginning
%     b = min( find(Alignment(1,:)~='-',1,'last'), find(Alignment(3,:)~='-',1,'last') ) ; % ignore gaps at end
%     % calculate Hamming distance
%     num_error(n) = hammingd( Alignment(1,:), a, b, Alignment(3,:), a, b ) ;
    num_error(n) = hammingd_trim_ends_gap(haplo(n), trueh(n)) ;
  end
end