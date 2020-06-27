%% MSA
clear; clc; close all

%% Check License
checkbio

%% Paths
path_out = './MATLAB/MSAga/data/out';
[status, msg, msgID] = mkdir(path_out);
filename_p53_multialign = fullfile(path_out,'p53_library.aln');

%% Loading FASTA Input
% This file should be in MATLAB path
p53 = fastaread('p53samples.txt')

%% Multialign
SeqsMultiAligned = multialign(p53)

%% Show
multialignwrite(filename_p53_multialign, SeqsMultiAligned)
S = multialignread(filename_p53_multialign)
seqalignviewer(S);

%% Fitness MSA
score = fitness_msa(S);
fprintf("Fitness Score = %.2f\n", score);
% Fitness Score = 7159.00