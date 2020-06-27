function [max_length] = calc_max_length(sequences_table)
num_sequences = size(sequences_table, 1);
max_length = 0;
for i = 1:num_sequences
    length = strlength(sequences_table{i,1}{1});
    if length > max_length
        max_length = length;
    end
end
end

