function [sum_score] = fitness_msa(multialign_read)
    sum_score = 0;
    for i = 2:size(multialign_read,2)                                
        % Compute pairwise scores
        for j = 1:(i-1)
            [curr_score, ~] = nwalign(multialign_read(i).Sequence, ...
                                      multialign_read(j).Sequence, ...
                                      'ScoringMatrix', 'BLOSUM62', ...
                                      'GapOpen', 1, ...
                                      'ExtendGap', 0.5);
            sum_score = sum_score + curr_score;
        end
    end
end