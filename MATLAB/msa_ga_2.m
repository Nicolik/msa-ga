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

%% Creating object
gen_alg = msaga(chromosomes, generations, min_num_gen, ....
                mutation_rate, crossover_prob, VERBOSE);
fprintf(gen_alg.toString());

%% Optimization Run
[align_cell, pop, best_chromosomes, stats] = ...
    gen_alg.run_ga(input_file, isFasta);

%% Convert To Struct
seq_struct = cell_to_struct(input_file, align_cell);

%% Show
multialignwrite(filename_p53_multialign, seq_struct)
S = multialignread(filename_p53_multialign);
seqalignviewer(S);

%% Fitness MSA
score = fitness_msa(S);
fprintf("Fitness Score = %.2f\n", score);
% Fitness Score = 7165.00
% Fitness Score = 7179.50
% Fitness Score = 7181.25

%% Optimization Results
fprintf("Alignment Result:\n");
disp(align_cell);

%% Plot Trends
id = 2;
plot_trends(stats, id, gen_alg.toString())
