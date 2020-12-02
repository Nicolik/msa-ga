classdef msaga
    %MSAGA Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        chromosomes
        generations
        min_generations
        mutation_rate
        crossover_prob
        VERBOSE
    end
    
    methods
        function obj = msaga(chromosomes, generations, min_generations, mutation_rate, crossover_prob, VERBOSE)
            %MSAGA Construct an instance of this class
            %   Detailed explanation goes here
            obj.chromosomes = chromosomes;
            obj.generations = generations;
            obj.min_generations = min_generations;
            obj.mutation_rate = mutation_rate;
            obj.crossover_prob = crossover_prob;
            obj.VERBOSE = VERBOSE;
        end
        
        function [str] = toString(obj)
            str = "";
            str = str + sprintf("+=========================================+\n");
            str = str + sprintf("+================ MSAGA ==================+\n");
            str = str + sprintf("+=========================================+\n");
            str = str + sprintf("+== Chromosomes           ==     %3d    ==+\n", obj.chromosomes);
            str = str + sprintf("+== Generations (Min)     ==     %3d    ==+\n", obj.min_generations);
            str = str + sprintf("+== Generations (Max)     ==     %3d    ==+\n",  obj.generations);
            str = str + sprintf("+== Crossover Probability ==    %.2f    ==+\n", obj.crossover_prob);
            str = str + sprintf("+== Mutation Rate         ==    %.2f    ==+\n", obj.mutation_rate);
            str = str + sprintf("+=========================================+\n");
        end    
        
        function [align_cell, pop, best_chromosomes, stats] = ...
                run_ga(obj, input_path, isFasta)
            % Read the input file sequences
            [seq_table] = prepare_input(input_path, isFasta);
            
            % Show the original chromosome
            fprintf("Input Sequence Table:\n");
            disp(seq_table);
            
            % Initialize the population
            pop = obj.initialize_pop(seq_table);       
            
            % Repeat for all generations or until a good solution is found
            best_val = false;
            best_chromosome = false;
            best_chromosomes = {};
            stats.best_values = [];
            stats.avg_values = [];
            stats.worst_values = [];
            stats.chromosomes = obj.chromosomes;
            gener_count = 1;
            new_pop = {};
            
            exec_time = 0;
            
            while gener_count < obj.generations
                
                tic
                
                % Evaluate all chromosomes
                evaluations = {};
                for i = 1:size(pop,1)
                    pop{i,1}.evaluation = msaga.eval_function(pop{i,1}.chromosome);
                    evaluations{end+1,1} = pop{i,1}.evaluation;
                end
                total_eval_time = toc;
                
                % Repeat for all chromosomes
                for chromosome_counter = 1:obj.chromosomes
                    % Get crossover probability
                    crossover_prob_ = rand();
                
                    % Crossover check
                    if crossover_prob_ < obj.crossover_prob
                        % Select two parents
                        [p1, p2] = msaga.select_parents(pop, evaluations, obj.VERBOSE);
                        
                        % Apply a crossover operation on p1 and p2
                        child = msaga.apply_crossover(p1, p2);
                        
                        % Apply a mutation on child
                        child.chromosome = obj.apply_mutation(child.chromosome);
                       
                        % Add Gaps (all chromosomes must be of same length)
                        child.chromosome = add_gaps(child.chromosome);
                        
                        % Remove columns composed of only gaps
                        child.chromosome = remove_useless_gaps(child.chromosome);
                        
                        % Add the child to the new population
                        new_pop{end+1,1} = child;
                        
                        % Display the child
                        if obj.VERBOSE
                            fprintf("\nChild:\n");
                            disp(child.chromosome);
                        end
                    end
                end
                
                % Get the best chromosome
                [best_chromosome, avg_val, best_val, worst_val] = ...
                    msaga.get_best_chromosome(new_pop);
                                
                % Add best chromosome to list
                best_chromosomes{end+1,1} = best_chromosome;
                
                % Stats
                stats.best_values(end+1,1) = best_val;
                stats.avg_values(end+1,1) = avg_val;
                stats.worst_values(end+1,1) = worst_val;
                
                % Break the execution if there are no relevant changes
                if obj.check_no_changes(stats.best_values)
                    fprintf("Stop the optimization. No progresses!\n");
                    break;
                end
                
                % Update population
                % Add to the population the new generated chromosomes and
                % remove the same number of the worst chromosomes from
                % the original population                
                [pop] = msaga.update_population(pop, new_pop);
                
                new_pop = {};                            
                end_time = toc;
                 % Print time and stats
                fprintf("[Gen %03d] Eval / Tot: %.2f / %.2f sec --  Worst / Avg / Best: %.1f / %.1f / %.1f \n", ...
                    gener_count, total_eval_time, end_time, worst_val, avg_val, best_val);
                exec_time = exec_time + end_time;
                gener_count = gener_count + 1;
            end
            fprintf("[Whole optimization] Elapsed Time: %.4f seconds\n", exec_time);
            [stats.best_value, stats.best_gen] = max(stats.best_values);
            fprintf("Best Result -- Gen %03d -- Val %.1f\n", ...
                stats.best_gen, stats.best_value);
            align_cell = best_chromosomes{stats.best_gen,1};
        end
        
        function [change] = check_no_changes(obj, best_values)
            num_best_chromosomes = size(best_values, 1);
            if num_best_chromosomes < obj.min_generations
                change = false;
                return;
            else
                percent = round(0.2 * num_best_chromosomes);
                if percent < 2
                    percent = 2;
                end
                last_best_values = best_values(end-percent:end,1);
                var_lbv = var(last_best_values);
                if var_lbv < 1.05
                    change = true;
                    return;
                end
            end
            change = false;
        end    
            
        function [child] = apply_mutation(obj, child)
            n = size(child,1);
                     
            
            % Select mutation method and apply
            if rand() < obj.mutation_rate
                if obj.VERBOSE
                    fprintf("\nChild before mutation:\n");
                    disp(child);
                end
                
                rand_value = rand();
                
                % Apply gaps removal
                if rand_value < 0.5
                    cell_i = randi([1,n]);
                    cell_j = randi([1,strlength(child{cell_i,1})]);
                    
                    child_i = char(child{cell_i,1});
                    child_i_j = child_i(cell_j);
                    while child_i_j ~= "-"
                        cell_i = randi([1,n]);
                        cell_j = randi([1,strlength(child{cell_i,1})]);
                        child_i = char(child{cell_i,1});
                        child_i_j = child_i(cell_j);
                    end
                    
                    % Get extending gap
                    gaps = get_interval_gaps(child, cell_i, cell_j);
                    start_gap = gaps{1,1};
                    end_gap = gaps{end,1};
                    if start_gap ~= end_gap && obj.VERBOSE
                        fprintf("Start gap = %d\tEnd gap = %d\n", start_gap,end_gap);                   
                    end
                        
                    % Remove gaps
                    start_end_child = char(child{cell_i,1});
                    start_child = start_end_child(1:start_gap);           
                    end_child = start_end_child(end_gap+1:end);
                    child{cell_i,1} = string(start_child) + string(end_child);                 
                    
                % Apply k condition (add gaps)
                else
                    cell_i = randi([1,n]);
                    cell_j = randi([1,strlength(child{cell_i,1})]);
                    k = randi(1, ceil(0.1 * calc_max_length(child)));
                    to_add = "";
                    
                    for i = 1:k
                        to_add = to_add + "-";
                    end
                    
                    % child[cell_i] = child[cell_i][:cell_j] + to_add + child[cell_i][cell_j:]
                    start_end_child = char(child{cell_i,1});
                    start_child = start_end_child(1:cell_j);           
                    end_child = start_end_child(cell_j+1:end);
                    child{cell_i} = string(start_child) + to_add + string(end_child);
                    
                end
                if obj.VERBOSE
                    fprintf("\nChild after mutation:\n");
                    disp(child);
                end
            end
            
        end
        
        function [pop] = initialize_pop(obj, seq_table)
            pop = {};
            for chromosome_counter = 1:obj.chromosomes
                
                % Create new chromosome
                new_chromosome = {};
                
                % Use nwalign to compute pairwise alignments
                % by the Needleman-Wunsch algorithm
                for i = 1:size(seq_table,1)
                    alignments = {};
                    
                    % Compute pairwise alignments
                    for j = 1:size(seq_table,1)
                        if i ~= j
                            [current_score, current_align] = nwalign(seq_table{i,1}{1}, seq_table{j,1}{1});
                            alignments{end+1,1} = string(current_align);
                        end
                    end
                    
                    % Randomly select an alignment
                    alignment = randsample(alignments, 1);
                    alignment = alignment{1}(1,:);
                    new_chromosome{end+1,1} = alignment;
                end
                
                % Add the chromosome generated and prints it
                new_chromosome = add_gaps(new_chromosome);
                pop_iter.chromosome = new_chromosome;
                pop_iter.evaluation = 0;
                pop{end+1,1} = pop_iter;
                if obj.VERBOSE
                    fprintf("\nChromosome %3d on %3d:\n", chromosome_counter, obj.chromosomes);
                    disp(new_chromosome);
                end
            end
        end
               
    end
    
    methods(Static)
        function [sum_score] = eval_function(chromosome)
            sum_score = 0;
            for i = 2:size(chromosome,1)                                
                % Compute pairwise scores
                for j = 1:(i-1)
                    [curr_score, ~] = nwalign(chromosome{i,1}, chromosome{j,1}, ...
                                              'ScoringMatrix', 'BLOSUM62', ...
                                              'GapOpen', 1, ...
                                              'ExtendGap', 0.5);
                    sum_score = sum_score + curr_score;
                end
            end
        end
        
        function [p1, p2] = select_parents(pop, evaluations, VERBOSE)
            % Normalize the fitness
            pop_sum = sum(cell2mat(evaluations));
            evaluations_aux = {};
            for i = 1:size(evaluations,1)
                evaluations_aux{i,1} = evaluations{i,1} / pop_sum;
            end
            
            % Build the mating pool
            pool = {};
            for i = 1:size(pop,1)
                prob = ceil(evaluations_aux{i,1} * 100);
                for j = 1:prob
                    pool{end+1,1} = pop{i,1};
                end
            end
            
            % Randomly select two parents
            p1 = randsample(pool,1);
            p1 = p1{1};
            p2 = randsample(pool,1);
            p2 = p2{1};
            
            % Show chromosomes
            if VERBOSE
                fprintf("\nParent 1:\n");
                disp(p1.chromosome);
                fprintf("\nParent 2:\n");
                disp(p2.chromosome);
            end
        end
        
        function [child_struct] = apply_crossover(p1, p2)
            n = numel(p1.chromosome);
            
            % Apply horizontal crossover
            rand_horizontal = randi([1,n-1]);
            child = {};
            for i = 1:n
                if i <= rand_horizontal
                    child{i,1} = p1.chromosome{i,1};
                else 
                    child{i,1} = p2.chromosome{i,1};
                end
            end
            child_struct.chromosome = child;
            child_struct.evaluation = 0;
        end
        
        function [best_chromosome, avg_val, best_val, worst_val] = get_best_chromosome(new_pop)
                avg_val = 0;
                best_val = 0;
                worst_val = Inf;
                best_chromosome = false;
                for i = 1:size(new_pop,1)
                    curr_val = msaga.eval_function(new_pop{i,1}.chromosome);
                    new_pop{i,1}.evaluation = curr_val;
                    avg_val = avg_val + curr_val;
                    if curr_val >= best_val
                        best_val = curr_val;
                        best_chromosome = new_pop{i,1}.chromosome;
                    end
                    if curr_val <= worst_val
                        worst_val = curr_val;
                    end
                end
                avg_val = avg_val / size(new_pop,1);
        end
        
        function [pop] = update_population(pop, new_pop)
                pop_evaluations = [];
                for i = 1:size(pop,1)
                    pop_evaluations = [pop_evaluations, pop{i,1}.evaluation];
                end
                [~, pop_indices] = sort(pop_evaluations, 'descend');
                new_pop_evaluations = []; 
                for i = 1:size(new_pop,1)
                    new_pop_evaluations = [new_pop_evaluations, new_pop{i,1}.evaluation];
                end
                [~, new_pop_indices] = sort(new_pop_evaluations, 'descend');
                
                for i = 1:size(new_pop,1)
                    index_pop = pop_indices( numel(pop_indices)-numel(new_pop_indices)+i );
                    index_new_pop = new_pop_indices (i);
                    pop{index_pop,1} = new_pop{index_new_pop,1}; 
                end 
        end
        
    end    
end

