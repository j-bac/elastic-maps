load('healthyTissues.mat');
[d3efn, col] = degap(d3efn, col);
stret1 = [0.1, 0, 0, 0, 0,];
bend1 = [0.1, 0.1, 0.01, 0.2, 0.5];
stret2 = stret1(1:4);
bend2 = [1, 1, 0.5, 0.1];
colors = ['k', 'm', 'g', 'b', 'm', 'g', 'b', 'm', 'g', 'b', 'm', 'g'];
shapes = ['s', 's', 's', 'o', 'o', 'o', 'd', 'd', 'd', '^', '^', '^'];
%all figures are written into subfolder tmp.
tim1 = oneTest(d3efn, [10, 10], stret1, bend1, {},{'classes', colN, 'markColour', colors, 'markShape', shapes});
%tim2 = oneTest(d3efn, [10, 10], stret2, bend2, {},{'classes', colN, 'markColour', colors, 'markShape', shapes});