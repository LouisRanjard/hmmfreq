function y = list_combi(n)
%% list the n choose k with repetition and order matters, number of states if n^k
% k: number of true haplotypes with known frequencies (hardcoded to 3)
% n: number of subsequences with observed frequencies
% returns the subsequence number of the (3) haplotypes

if nargin==0
    n=4;
end

x = 1:n; % the elements you want to choose from
i = 1:length(x); % indices into the vector x
[i1,i2,i3] = ndgrid(i); % all possible combinations of indices into x
y = x([i3(:) i2(:) i1(:)]);


