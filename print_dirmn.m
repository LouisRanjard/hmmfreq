function [] = print_dirmn(obsc,regionc)
% print dirichlet multinomial probability for different combinations of the
% 3 observed counts given local regional background count

fprintf(1,'obs:\t');
for a=1:length(obsc)
    fprintf(1,'%d\t',obsc(a)) ;
end
for a=1:length(obsc)
    fprintf(1,'%.2f%%\t',100*obsc(a)/sum(obsc)) ;
end
fprintf(1,'\n');

fprintf(1,'loc:\t');
for a=1:length(regionc)
    fprintf(1,'%d\t',regionc(a)) ;
end
for a=1:length(regionc)
    fprintf(1,'%.2f%%\t',100*regionc(a)/sum(regionc)) ;
end
fprintf(1,'\n');

y = list_combi(3) ;
exp_pct = [62.5 25 12.5] ;
%exp_count = zeros(1,3) ;
for m=1:size(y,1)
    fprintf(1,'%2d| %d %d %d',m,y(m,:)) ;
    counts = zeros(1,size(obsc,2)) ;
    model = zeros(1,size(obsc,2)) ;
    unik = unique(y(m,:)) ;
    for m2 = unik
        for m1 = 1:size(obsc,2)
            counts(m2) = counts(m2) + regionc(m1)*(y(m,m1)==m2) ;
            model(m2) = model(m2) + exp_pct(m1)*(y(m,m1)==m2) ;
            %exp_count(m2) = exp_count(m2) + exp_pct(m1)*(y(m,m1)==m2) ;
        end
    end
    counts(counts==0) = round(0.05*sum(counts)) ; % allow 5% error in counts
    fprintf(1,'\t%.2f',dirmnpdfln(obsc,counts));
    %fprintf(1,'\t%.2f',dirmnpdfln(obsc,exp_pct));
    fprintf(1,'\t\t');
    fprintf(1,'%.1f\t',model);
    fprintf(1,'\n');
end

end