function [seq_table] = prepare_input(input_path, isFasta)
if isFasta
    seq_table = fastaread_as_table(input_path);
else
    seq_table = readtable(input_path, 'ReadVariableNames', 0);
end
seq_table.Properties.VariableNames = "Sequences";
num_sequences = size(seq_table, 1);
for i = 1:num_sequences
    seq_table{i,1}{1} = string(seq_table{i,1}{1});
end
[seq_table] = add_gaps(seq_table);
end
