We first ran both algorithms with default values (implementations from clusterMaker2).

Then we tried MCODE with cutoff = 0.1 instead of 0.2, which reduced number of nodes clustered.
On the other hand, increasing it to 0.4 also increased number of nodes.
We also reduced number of iterations for MCL from 16 to 4, partially because the original value meant quite long computations, compared to MCODE especially. It reduced number of clusters :)
We then added SCPS edge weight modification on top of reduced number of iterations. Nothing happened though, maybe there were no weights to begin with.

Also, we ran MCODE using some other package, that generated a whole analysis included in separate .txt file.

The results of MCODE and MCL are incomparable, MCL clusters almost everyting, while MCODE dropped most of the nodes. We tried making these more similar, but they were simply too far from each other.