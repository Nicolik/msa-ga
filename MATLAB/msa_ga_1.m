%% MSAGA - Example 1
clear; clc; close all

%% Check License
checkbio

%% Configuration
chromosomes = 32;
generations = 200;
min_num_gen = 80;
mutation_rate = 0.1;
crossover_prob = 0.5;
isFasta = false;
VERBOSE = false;

%% Creating object
input_file = "data_ex_1.txt";
gen_alg = msaga(chromosomes, generations, min_num_gen, ....
                mutation_rate, crossover_prob, VERBOSE);
fprintf(gen_alg.toString());

%% Optimization Run
[align_cell, pop, best_chromosomes, stats] = ...
    gen_alg.run_ga(input_file, isFasta);

%% Optimization Results
fprintf("Alignment Result:\n");
disp(align_cell);

%% Plot Trends
id = 1;
plot_trends(stats, id, gen_alg.toString())