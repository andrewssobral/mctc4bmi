# mctc4bmi
Matrix and Tensor Completion for Background Model Initialization

### Usage example
Please refer to [main.m](https://github.com/andrewssobral/mctc4bmi/blob/master/main.m) and [demo.m](https://github.com/andrewssobral/mctc4bmi/blob/master/demo.m) files.

### Citation
If you use this code for your publications, please cite it as:
```
@article{sobral_sbmi_prl_2016,
	title   = {Matrix and Tensor Completion Algorithms for Background Model Initialization: A Comparative Evaluation},
	author  = {Sobral, A. and Zahzah, E.},
	journal = {Pattern Recognition Letters},
	doi     = {10.1016/j.patrec.2016.12.019},
	issn    = {01678655},
	year    = {2016}
}
```

### Methods
List of the algorithms:

----------------------------------------------
|<sub>Type </sub>|<sub>Algorithm name <br/>(click to see the source code) </sub>|<sub> Author(s) </sub>|
|----------------|--------------------------------------------------------------|----------------------|
|<sub>**Matrix Completion**</sub>| | |
|                | [GROUSE](https://github.com/andrewssobral/mctc4bmi/blob/master/algs_mc/GROUSE) | [Balzano and Wright (2013)](http://ieeexplore.ieee.org/document/6713992/) |
|                | [IALM](https://github.com/andrewssobral/mctc4bmi/blob/master/algs_mc/IALM) | [Lin et al.(2010)](https://arxiv.org/abs/1009.5055) |
|                | [LMaFit](https://github.com/andrewssobral/mctc4bmi/blob/master/algs_mc/LMaFit) | [Wen et al.(2012)](http://link.springer.com/article/10.1007/s12532-012-0044-1) |
|                | [LRGeomCG](https://github.com/andrewssobral/mctc4bmi/blob/master/algs_mc/LRGeomCG) | [Vandereycken (2013)](https://arxiv.org/abs/1209.3834) |
|                | [MC-NMF](https://github.com/andrewssobral/mctc4bmi/blob/master/algs_mc/MC-NMF) | [Xu et al.(2012)](https://arxiv.org/abs/1103.1168) |
|                | [OptSpace](https://github.com/andrewssobral/mctc4bmi/blob/master/algs_mc/OptSpace) | [Keshavan et al.(2010)](https://arxiv.org/abs/0906.2027) |
|                | [OR1MP](https://github.com/andrewssobral/mctc4bmi/blob/master/algs_mc/OR1MP) | [Wang et al.(2015)](https://arxiv.org/abs/1404.1377) |
|                | [RMAMR](https://github.com/andrewssobral/mctc4bmi/blob/master/algs_mc/RMAMR) | [Ye et al.(2015)](http://ieeexplore.ieee.org/abstract/document/7014298/) |
|                | [ScGrassMC](https://github.com/andrewssobral/mctc4bmi/blob/master/algs_mc/ScGrassMC) | [Ngo and Saad (2012)](http://dl.acm.org/citation.cfm?id=2999292) |
|                | [SVP](https://github.com/andrewssobral/mctc4bmi/blob/master/algs_mc/SVP) | [Jain et al.(2010)](https://arxiv.org/abs/0909.5457) |
|<sub>**Tensor Completion**</sub>| | |
|                | [TMac](https://github.com/andrewssobral/mctc4bmi/tree/master/algs_tc/TMac) | [Xu et al. (2015)](https://arxiv.org/abs/1312.1254) |
