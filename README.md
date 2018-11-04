# Building Optimized-Genetic-Codes with Neural Networks and Evolving HIV strains with Genetic Algorithms
Our standard genetic code is built so mutations, that naturally occur, won't take down our whole system, i.e. we're not too big to fail. A long standing problem has been to develop other genetic codes that are more resistent to the effects of mutation.  The files on this page contain an AI,using Ver Der Bout's (1989) tweaked Hopfield Network algorithm, to generate genetic codes that functions well under pressures of mutation.  Further, the AI is flexible allowing you to input parameters you'd like to optimize your code for, e.g. Polarity, Hydropathy, Volume, Iso-Electricity. This flexibility separates this method of generating genetic codes from others.

A full description of the work and results are given in the preprint below-

"Phylogeny Recapitulates Learning: Self-Optimization of Genetic Code" (Attie,Sulkow,Di,and Qiu,2018)
https://www.biorxiv.org/content/early/2018/02/06/260877

The full set of files associated with this paper can be found at https://github.com/weigangq/code-by-tsp.  

aaPaths2020- The file contains the AI algorithm that's used to find paths through the set of 20 amino acids. These paths can be optimized for any combination of dimensions that quantitatively describe Amino Acids. For this program we have restricted the 4 dimensions- Polarity, Hydropathy, Volume and Iso-Electricity.

geneticAlgorithmsforHIV- Using genetic algorithms we seek to evolve sets of HIV alleles, represented by binary sequences, so that minimum distance between any pair of alleles.This question is relevant to immunologists looking for similarities between HIV strains that aren't related by recombinations. It's also relevant to computer science folks looking for to maximize minimum distance between codewords with a prescribed length. 

pairplots2-  PairPlots2 gives you a method for comparing the distribution of errors associated with a genetic code. 

2020Plots-This 2D plot, created with GGPLOT2, shows the distributions of errors associated to genetic codes created with paths through 20 amino acids.

2121Plots- This 2D plot, created with GGPLOT2, shows the distributions of errors associated to genetic codes created with paths through 20 amino acids.

pcaPlots- This file allows you to graph paths through the amino acids using coordinates (polarity, hydropathy, iso-electricity, volume) with amino acids and PCA coordinates.   

pcaPlots.pdf- This 2D plot, created with GGPLOT2,  shows paths through amino acids using PCA coordinates.  One can visual compare random path through 21 amino acids, with a path induced by the standard genetic code, and an optimized code generated with our AI.  
