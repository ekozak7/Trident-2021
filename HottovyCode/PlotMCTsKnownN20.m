% Script to print MCTS with CIs 

Loops = [100 200 400 800 1600 3200 6400 12800];
Means = [19.21 11.29 6.54 3.87 2.03 1 0.5 0.19];
Vars = [579.783 171.569 69.956 26.814 9.501 3.315 1.497 0.482];

LoopConf = [Loops Loops(end:-1:1)];
MeanConf = [Means+sqrt(Vars) Means(end:-1:1)-sqrt(Vars(end:-1:1))];

figure
p = fill(LoopConf, MeanConf,'blue');
p.FaceColor=[0.8 0.8 1];
p.EdgeColor = 'none';

hold on
plot(Loops, Means, 'ko')
hold off