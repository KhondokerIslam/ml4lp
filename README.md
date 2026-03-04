<!-- <div align="center"> -->
<!-- omit in toc -->
# Evaluating SoTa PLMs on Long Sequence Generation

<a href="https://pytorch.org/get-started/locally/"><img alt="PyTorch" src="https://img.shields.io/badge/PyTorch-ee4c2c?logo=pytorch&logoColor=white"></a>
<a href="https://pytorchlightning.ai/"><img alt="Lightning" src="https://img.shields.io/badge/-Lightning-792ee5?logo=pytorchlightning&logoColor=white"></a>
<a href="https://hydra.cc/"><img alt="Config: Hydra" src="https://img.shields.io/badge/Config-Hydra-89b8cd"></a>
<a href="https://github.com/ashleve/lightning-hydra-template"><img alt="Template" src="https://img.shields.io/badge/-Lightning--Hydra--Template-017F2F?style=flat&logo=github&labelColor=gray"></a>

## Overview
This repository systematically evaluates [CFP-GEN (Combinatorial Functional Protein Generation via Diffusion Language Model)](https://arxiv.org/pdf/2505.22869) on **long protein sequence generation (>1,000 amino acids)**.

While CFP-GEN achieves state-of-the-art results on standard protein generation benchmarks, its ability to generalize to long multi-domain proteins remains unexplored. This project investigates:

* Does multi-constraint conditioning help long-sequence generation?
* Does performance degrade as sequence length increases?
* How does CFP-GEN compare to its unconditional backbone (DPLM)?

This work is built directly on top of the official CFP-GEN implementation.

* 📄 CFP-GEN Paper: [https://arxiv.org/pdf/2505.22869](https://arxiv.org/pdf/2505.22869)
* 💻 CFP-GEN Repository: [https://github.com/yinjunbo/cfpgen](https://github.com/yinjunbo/cfpgen)
* 📝 My term paper: [./asset/term_paper.pdf](./asset/term_paper.pdf)


## Abstract 

 Protein language models (PLMs) have achieved remarkable progress in discriminative and generative protein modeling, especially under biologically meaningful constraints such as InterPro (IPR) domains and Gene Ontology (GO) terms. Multi-constraint diffusion-based frameworks such as CFP-GEN demonstrate state-of-the-art performance on standard benchmarks. However, their ability to generalize to **long protein sequence generation (>1,000 residues)** remains largely unexplored, despite the biological relevance of large multi-domain proteins. In this work, I systematically evaluate CFP-GEN on long-sequence conditional generation using a curated Swiss-Prot dataset (1,000–2,000 residues). I compare it against its unconditional backbone (DPLM) under controlled IPR constraints.

 Evaluation using Mean Reciprocal Rank (MRR), and Maximum Mean Discrepancy (MMD), reveals that although CFP-GEN maintains comparable ranking performance, it does not outperform—and sometimes underperforms—the unconditional baseline in distributional alignment for long sequences. Performance degradation increases with sequence length for both models. These findings suggest that incorporating biological constraints alone does not resolve long-context generation challenges in PLMs.

## CFP-Gen 🌟
Please refer to CFP-GEN [paper](https://arxiv.org/pdf/2505.22869) and their respective [repository](https://github.com/yinjunbo/cfpgen) for better more context of this task. Also refer to my [term paper](./asset/term_paper.pdf) for understanding my work.


## Installation

```bash
# clone project
git clone --recursive https://github.com/KhondokerIslam/ml4lp.git

# create conda virtual environment
env_name=cfpgen

conda create -n ${env_name} python=3.9 pip
conda activate ${env_name}

# automatically install everything else
bash install.sh

# alternative
pip install -r  requirement.txt
```


## Datasets
**CFP-GEN General Dataset (Paper Dataset)**: The processed dataset ```cfpgen_general_dataset``` can be downloaded from [Google Drive](https://drive.google.com/file/d/1bRtil483NBOuazPSVO7gCpM-K7rNj1Z9/view?usp=sharing) and placed in the directory: ```dataset/general_dataset_go_ipr/cfpgen_general_dataset/```. It contains 103,939 proteins annotated with 375 GO terms and 1,154 IPR domains.

**GO/IPR Mapping**: CFP-GEN only supports generation following ground-truth labels from natural proteins (e.g., from SwissProt).
The [GO/IPR mapping info](https://drive.google.com/drive/folders/1Z6Zmjy1h41rk_Lu89itHeWdep-S9-Zjy?usp=sharing) can be download here. Plase place them on ```dataset/ipr_mapping.json```.

**Long Sequence Dataset (This Work)**: We evaluate on proteins between **1,000 and 2,000 amino acids**. You can download long Swiss-Prot sequences from [UniProt](https://www.uniprot.org/uniprotkb?query=*&facets=reviewed%3Atrue%2Clength%3A%5B801+TO+*%5D). Place the downloaded ```.fasta``` and ```.tsv``` files inside ```dataset/```. For reproducibility, the dataset used in this study is already included: ```ataset/uniprotkb_AND_reviewed_true_AND_model_o_2026_02_18.fasta```, ```dataset/uniprotkb_AND_reviewed_true_AND_model_o_2026_02_18.tsv```

**Generate Final Evaluation Set**: Run `generate_test_suit.py`. This will generate ```data-bin/uniprotKB/cfpgen_general_dataset/experiment.pkl``` required to run this experiment.



### Notes:

- ```cfpgen-650m```: Support conditioning on GO terms, IPR domains and sequence motifs (e.g., 10-30 residue fragments) defined by our **general protein dataset**. This model can be readily used for _Functional Protein Generation_.

- ```dplm-650m```: This is the base pretrained model from DPLM, required to be placed under ```cfpgen/pretrained/```.


## Generation with _CFP-Gen_
### Functional Protein Generation

Users could modify necessary parameters (e.g.,```ckpt_path=<path_to_cfpgen-650m>```) in the config file:
```bash
configs/test_cfpgen.yaml
```
and then run the following command to start generation:

```bash
python cfp_generate.py
```
The results will be saved in `./generation-results`.


# Evaluation
We only evaluate on sequence level (i.e., distributional statistics) and in particularly on IPR domains.

### Distribution Evaluation
The following command computes Maximum Mean Discrepancy (MMD) and Mean Reciprocal Rank (MRR) between the generated and real sequences:
```bash
python eval_mmd.py <ipr> <fasta_filename> <gt_data>
```
Here, ```<fasta_filename>```is the output FASTA file obtained by the previous generation commands. ```<gt_data>``` refers to the ground-truth data file (e.g., ```data-bin/uniprotKB/cfpgen_general_dataset/test.pkl```)

# Citation

```
To be added soon.
```
