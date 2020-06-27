function [seq_table] = fastaread_as_table(input_path)
seq_struct = fastaread(input_path);
seq_table = struct2table(seq_struct);
seq_table = seq_table(:,2);
end

