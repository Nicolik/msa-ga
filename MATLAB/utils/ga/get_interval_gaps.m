function [gaps] = get_interval_gaps(lines, cell_i, cell_j)
aux_cell_j = cell_j;
gaps = {};
symbol_found = false;

% Find gaps after cell_j (inclusive)
while ~symbol_found
    if cell_j > strlength(lines{cell_i,1})
        break;
    else
        lines_i = char(lines{cell_i,1});
        lines_i_j = lines_i(cell_j);
        if lines_i_j == "-"
            gaps{end+1,1} = cell_j;
            cell_j = cell_j + 1;
        else 
            symbol_found = true;
        end 
    end
end

% Find gaps before cell_j (exclusive)
aux_cell_j = aux_cell_j - 1;
while ~symbol_found
    if aux_cell_j < 1
        break;
    else
        lines_i = char(lines{cell_i,1});
        lines_i_j = lines_i(aux_cell_j);
        if lines_i_j == "-"
            gaps = [ { aux_cell_j }; gaps ];
            aux_cell_j = aux_cell_j - 1;
        else
            symbol_found = true;
        end        
    end
end

end

