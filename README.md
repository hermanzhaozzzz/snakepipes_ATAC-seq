# snakepipes_ATAC-seq
A standard ATAC-seq snakemake pipeline
## env
```shell
conda env create -f conda_env.yml
```
## run
```shell
# Jupyterlab运行step1生成json文件或者命令行运行笔记本
runipy step.01.GetFileName.ipynb
# 测试运行
snakemake -pr -j 10 -s step.02.Snakefile.py -n
# 实际运行
snakemake -pr -j 10 -s step.02.Snakefile.py
```
