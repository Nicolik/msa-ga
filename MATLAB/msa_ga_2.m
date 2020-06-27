%% MSA
clear; clc; close all

%% Check License
checkbio

%% Paths
path_out = './MATLAB/MSAga/data/out';
[status, msg, msgID] = mkdir(path_out);
filename_p53_multialign = fullfile(path_out,'p53_ga.aln');

%% Loading FASTA Input
% This file should be in MATLAB path
input_file = 'p53samples.txt';

%% Multialign with GA
chromosomes = 16;
generations = 200;
min_num_gen = 100;
mutation_rate = 0.05;
crossover_prob = 0.5;
isFasta = true;
VERBOSE = false;

gen_alg = msaga(chromosomes, generations, min_num_gen, ....
                mutation_rate, crossover_prob, VERBOSE);
[align_cell, pop, best_chromosomes, best_values] = ...
    gen_alg.run_ga(input_file, isFasta);

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