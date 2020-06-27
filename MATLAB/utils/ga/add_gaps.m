function [sequences_table] = add_gaps(sequences_table)
num_sequences = size(sequences_table, 1);
max_length = calc_max_length(sequences_table);
for i = 1:num_sequences
    sequence = sequences_table{i,1}{1};
    length = strlength(sequence);
    if length < max_length
        diff_length = max_length - length;
        for j = 1:diff_length
            sequence = sequence + "-";    
        end
        if istable(sequences_table)
            sequences_table{i,1} = {sequence};
        else
            sequences_table{i,1} = sequence;
        end
    end
end
end

