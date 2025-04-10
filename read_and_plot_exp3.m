%% Add Paths
addpath(genpath('./Utilities/'))
%% Read data and plot
cond_Pade = csvread('./data/exp3/cond_Pade.csv');
cond_Pade = cond_Pade(:, 1:100);
cond_Taylor = csvread('./data/exp3/cond_Taylor.csv');
cond_Taylor = cond_Taylor(:, 1:100);
m_Pade = csvread('./data/exp3/m_Pade.csv');
m_Pade = m_Pade(:, 1:100);
m_Taylor = csvread('./data/exp3/m_Taylor.csv');
m_Taylor = m_Taylor(:, 1:100);
succ_Pade = csvread('./data/exp3/succ_Pade.csv');
succ_Pade = succ_Pade(:, 1:100);
succ_Taylor = csvread('./data/exp3/succ_Taylor.csv');
succ_Taylor = succ_Taylor(:, 1:100);
data_plot(100, m_Taylor, m_Pade, cond_Taylor, cond_Pade, succ_Taylor, succ_Pade, './fig/fig3.eps')
%% Remove Path
rmpath(genpath('./Utilities/'))