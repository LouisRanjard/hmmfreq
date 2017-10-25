function [dhamming] = hammingd(s1,a1,b1,s2,a2,b2)
% return number of differences between two input sequences
% NB: function considers two 'N' as different
    if (isempty(s1) || isempty(s2))
        dhamming = b1-a1+1 ;
    else
        dhamming = sum( (s1(a1:b1) - s2(a2:b2))~=0 ) ; % positions that are different
        N2 = (s1(a1:b1)=='N') + (s2(a2:b2)=='N') ;     % positions where each sequence is 'N'
        dhamming = dhamming + sum(N2==2) ;
    end
end