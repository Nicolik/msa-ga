%% MSA
clear; clc; close all

%% Check License
checkbio

%% Paths
path_out = './MATLAB/data/out';
path_out_full = fullfile(pwd,path_out);
[status, msg, msgID] = mkdir(path_out);
filename_p53_multialign = fullfile(path_out,'p53_ga.aln');

%% Loading FASTA Input
% This file should be in MATLAB path
input_file = 'p53samples.txt';

%% Multialign with GA
chromosomes = 16;
generations = 200;
min_num_gen = 80;
mutation_rate = 0.2;
crossover_prob = 0.5;
isFasta = true;
VERBOSE = false;

gen_alg = msaga(chromosomes, generations, min_num_gen, ....
                mutation_rate, crossover_prob, VERBOSE);
[align_cell, pop, best_chromosomes, best_values, ...
    avg_values, worst_values] = gen_alg.run_ga(input_file, isFasta);

%% Convert To Struct
seq_struct = cell_to_struct(input_file, align_cell);

%% Show
multialignwrite(filename_p53_multialign, seq_struct)
S = multialignread(filename_p53_multialign)
seqalignviewer(S);

%% Fitness MSA
score = fitness_msa(S);
fprintf("Fitness Score = %.2f\n", score);
% Fitness Score = 7165.00
% Fitness Score = 7179.50

%% Optimization Results
best_value = max(best_values);
worst_value = min(worst_values);
std_avg = std(avg_values);
fprintf("Best value = %.1f\n", best_value);
fprintf("Alignment Result:\n");
disp(align_cell);

%% Plot Trends
figure('units','normalized','outerposition',[0 0 1 1])
plot(1:size(best_values,1), best_values,  '-o', ...
     1:size(best_values,1), avg_values,   ':*', ...
     1:size(best_values,1), worst_values, ':d', ...
     'LineWidth', 1.5)
ylim([worst_value - 2 * std_avg, best_value + 2 * std_avg])
title("Optimization Process Trends")
xlabel("Generation")
ylabel("Fitness")
legend("Best Values", "Average Values", "Worst Values")
saveas(gcf,'./MATLAB/images/trends_2.png')