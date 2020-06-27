function [seq_struct] = cell_to_struct(input_file, seq_cells)
seq_struct = fastaread(input_file);

num_sequences = size(seq_struct,1);
for i = 1:num_sequences
    seq_struct(i).Sequence=seq_cells{i,1}{1};
end
end

