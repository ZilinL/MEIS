# MEIS
This repository is the implementation of "Manifold embedded instance selection to suppress negative transfer in motorimagery-based brain–computer interface".

The MEIS algorithm operates in two ways: converting raw EEG matrices into manifold embedded vectors that maintain sample discriminability, and designing an evaluator to assess the transferability of samples and filter out
negative transfer samples from the source domain.
The flowchart of the proposed MEIS method is shown below.
![image](https://github.com/ZilinL/MEIS/assets/10232596/0be405b0-5bb3-4e70-98b5-72de19abe2bb)



## Datasets
1. [BCI competition IV-1](https://www.bbci.de/competition/iv/)
2. [BCI competition IV-2a](https://www.bbci.de/competition/iv/)
3. [SHU](https://figshare.com/articles/software/shu_dataset/19228725/1)
4. [OpenBMI](http://gigadb.org/dataset/view/id/100542)

## Running the code
The code works on Matlab 2022b software under Windows 11 system, and other versions of Matlab software should be able to run, but it has not been tested.

**DEMO_MEIS_MI1_SD.m:** A demo program that performs single-source domain experiments on the BCICIV-1 dataset.

**MEIS.m:** MEIS method code implementation.

**sample_reweight_kmm.m:** Incorporate KMM method for weight allocation to samples.

## Results
Experimental results and comparison methods are shown in the following table.
![image](https://github.com/ZilinL/MEIS/assets/10232596/3fcaabaf-6785-4edd-ba05-276093ac978a)

## Citation
If our work has been helpful to your research, please consider citing it.
```
@article{liang2024manifold,
  title={Manifold embedded instance selection to suppress negative transfer in motor imagery-based brain--computer interface},
  author={Liang, Zilin and Zheng, Zheng and Chen, Weihai and Pei, Zhongcai and Wang, Jianhua and Chen, Jianer},
  journal={Biomedical Signal Processing and Control},
  volume={88},
  pages={105556},
  year={2024},
  publisher={Elsevier}
}
```

## Reference
```
@Article{zhang2020manifold,
  title={Manifold Embedded Knowledge Transfer for Brain-Computer Interfaces},
  author={Zhang, Wen and Wu, Dongrui},
  journal={IEEE Trans. on Neural Systems and Rehabilitation Engineering},
  volume={28},
  number={5},
  pages={1117--1127},
  year={2020}
}
```

```
@article{huang2006correcting,
  title={Correcting sample selection bias by unlabeled data},
  author={Huang, Jiayuan and Gretton, Arthur and Borgwardt, Karsten and Sch{\"o}lkopf, Bernhard and Smola, Alex},
  journal={Advances in neural information processing systems},
  volume={19},
  year={2006}
}
```
