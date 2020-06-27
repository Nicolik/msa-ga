%% MSAGA - Example 1
clear; clc; close all

%% Check License
checkbio

%% Configuration
chromosomes = 16;
generations = 200;
min_num_gen = 100;
mutation_rate = 0.05;
crossover_prob = 0.5;
isFasta = false;
VERBOSE = false;

%% Creating object
input_file = "data_ex_1.txt";
gen_alg = msaga(chromosomes, generations, min_num_gen, ....
                mutation_rate, crossover_prob, VERBOSE);
            
%% Optimization Run
[align_cell, pop, best_chromosomes, best_values] = ...
    gen_alg.run_ga(input_file, isFasta);

%% Optimization Results
best_value = max(best_values);
worst_value = min(best_values);
fprintf("Best value = %d\n", best_value);
fprintf("Best chromosome\n");
disp(best_chromosomes{end,1});
fprintf("Alignment Result:\n");
disp(align_cell);
plot(best_values, 'o'), ylim([worst_value - 20, best_value + 20]), xlabel("Generation"), ylabel("Fitness")