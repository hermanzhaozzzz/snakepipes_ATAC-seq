# snakepipes_ATAC-seq
## about author

> author: [赵华男 | ZHAO Hua-nan](https://scholar.google.com/citations?user=ojSVoWQAAAAJ&hl=en)
>
> email: hermanzhaozzzz@gmail.com
>
> [Zhihu](https://www.zhihu.com/people/hymanzhaozzzz) | [BLOG](http://zhaohuanan.cc)
## doc
A standard ATAC-seq snakemake pipeline
## env
```shell
conda env create -f conda_env.yml
```
## run
```shell
# run Jupyter notebook to abtain the config
# run this cmd
# or
# open notebook and run all cells
runipy step.01.GetFileName.ipynb
# 测试运行
snakemake -pr -j 10 -s step.02.Snakefile.smk.py -n
# 实际运行
snakemake -pr -j 10 -s step.02.Snakefile.smk.py
```
## project structure