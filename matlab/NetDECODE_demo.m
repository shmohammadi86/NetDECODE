clear
addpath(genpath('utils'));
input_path = '../input/PBMC';


%% Import expression matrix
expression = mmread(fullfile(input_path, 'expression.mm'));
gene_names = readtable(fullfile(input_path, 'gene_names.txt'), 'ReadVariableNames', false);


%% Run NetDECODE
tic; G = NetDECODE(expression); toc
tic; G_pruned = sparse(pruneNetwork(G)); toc

deg = sum(G_pruned > 0);
ksdensity(deg)

A = expression > 0;
activity_scores = full(G_pruned'*A);

figure; plot(smooth(activity_scores(strcmp(gene_names.Var1, 'CD19'), :)))
figure; plot(smooth(activity_scores(strcmp(gene_names.Var1, 'CD3G'), :)))
figure; plot(smooth(activity_scores(strcmp(gene_names.Var1, 'MNDA'), :)))
