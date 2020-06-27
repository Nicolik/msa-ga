function [sequences_table] = remove_useless_gaps(sequences_table)
    to_rm = {};
    only_gaps = true;
    max_length = calc_max_length(sequences_table);
    
    % Find columns with gaps only
    for i = 1:max_length
        for j = 1:size(sequences_table,1)
            chromosome_j = char(sequences_table{j,1});
            chromosome_j_i = chromosome_j(i);
            if chromosome_j_i ~= " " && chromosome_j_i ~= "-"
                only_gaps = false;
            end
        end
        
        if only_gaps
            to_rm{end+1,1} = i;
        else
            only_gaps = true;
        end
    end
    
    % Remove gap-only columns
    to_rm = flip(to_rm);
    for i = 1:size(sequences_table,1)
        line = [];
        for j = 1:max_length
            chromosome_i = char(sequences_table{i,1});
            chromosome_i_j = chromosome_i(j);
            line = [line, chromosome_i_j];
        end
        
        for k = 1:numel(to_rm)
            to_rm_el = to_rm{k,1};
            line(to_rm_el) = [];
        end
        
        sequences_table{i,1} = string(line);
    end
end

