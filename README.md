[![Review Assignment Due Date](https://classroom.github.com/assets/deadline-readme-button-24ddc0f5d75046c5622901739e7c5dd533143b0c8e959d652212380cedb1ea36.svg)](https://classroom.github.com/a/Yyz_ETFO)
# ECE 208 Homework 6: Simulating Sequences Under GTR
In this homework assignment, we will study one of the robust models of DNA evolution: the [General Time Reversible (GTR) model](https://en.wikipedia.org/wiki/Models_of_DNA_evolution#GTR_model_(Tavar%C3%A9_1986)).<sup>1</sup> In this assignment, we will ignore insertions and deletions (indels), and we will assume that sequences only differ due to substitutions. We will assume that the evolutionary processes of the positions of a sequences are **identical and independently distributed (i.i.d.)**.

## Problem
The goal of this problem is to simulate the evolution of a sequence down a phylogenetic tree under the GTR model. Given a phylogenetic tree ***T*** with leafset ***L***, a root sequence ***s***, stationary probabilities **(*π*<sub>A</sub>, *π*<sub>C</sub>, *π*<sub>G</sub>, *π*<sub>T</sub>)**, and transition rates ***α* = *r*(A→G) = *r*(G→A)**, ***β* = *r*(A→C) = *r*(C→A)**, ***γ* = *r*(A→T) = *r*(T→A)**, ***δ* = *r*(G→C) = *r*(C→G)**, ***ε* = *r*(G→T) = *r*(T→G)**, and ***η* = *r*(C→T) = *r*(T→C)**, you will simulate the evolution of *s* down *T* to produce a sequence for each leaf in *L*.

For the root sequence *s*, you will need to be able to handle two scenarios:
1. The root sequence *s* may be given to you (e.g. [example/root_seq.fas](example/root_seq.fas))
2. You may be asked to generate a random sequence of length *k* to use as the root sequence *s*. In your report, you will describe how you simulate a sequence.

### Input
[generate_seqs.py](generate_seqs.py) will take as input the following items:
* `-t`: A phylogenetic tree (either binary or multifurcating) in the [Newick format](https://en.wikipedia.org/wiki/Newick_format)
* `-p`: A file containing the GTR parameters (see below)
* `-r`: Either the path to a [FASTA file](https://en.wikipedia.org/wiki/FASTA_format) containing the root sequence *s*, or an integer *k* representing the length of the root sequence *s* you will need to generate
* `-o`: The output file in the [FASTA format](https://en.wikipedia.org/wiki/FASTA_format) (default: standard output)

One function you will be modifying is the `evolve(tree, root_seq, gtr_probs, gtr_rates)` function, which takes as input the following parameters:
* `tree`: A rooted [TreeSwift](https://github.com/niemasd/TreeSwift) `Tree` object (either binary or multifurcating)
* `root_seq`: The root sequence *s*
* `gtr_probs`: A list of floats in the format `[prob_A, prob_C, prob_G, prob_T]` representing the stationary probabilities (sums to 1)
* `gtr_rates`: A list of floats in the format `[rate_CT, rate_AT, rate_GT, rate_AC, rate_CG, rate_AG]` representing the transition rates

The other function you will be modifying is the `random_seq(k, gtr_probs)` function, which takes as input the following parameters:
* `k`: An integer denoting the length of the random sequence to simulate
* `gtr_probs`: A list of floats in the format `[prob_A, prob_C, prob_G, prob_T]` representing the stationary probabilities (sums to 1)

#### GTR Parameter File
The GTR parameter file will be in the following format (an example can be found in [example/gtr_params.txt](example/gtr_params.txt):

```
prob_A prob_C prob_G prob_T
rate_CT rate_AT rate_GT rate_AC rate_CG rate_AG
```

**Note:** Recall that units of GTR rate matrix and units of branch length that you are scoring should be the same, or else the calculations will be wrong. One important note about the input GTR parameter file is that rates in this file **are NOT necessarily normalized to be in substitution units**. Thus, **normalization is required** and is left to your code. In the write-up, you will also need to describe how you normalize the GTR transition rates.

### Output
[generate_seqs.py](generate_seqs.py) will output a [FASTA file](https://en.wikipedia.org/wiki/FASTA_format) containing a simulated sequence for each leaf in the phylogenetic tree.

The `evolve(tree, root_seq, gtr_probs, gtr_rates)` function will return a dictionary where keys are the labels of the leaves of `Tree` (which are strings) and values are the corresponding sequences (which are also strings).

The `random_seq(k, gtr_probs)` function will return a string of length `k`.

### Coding
* To perform operations on matrices, including matrix exponential and matrix logarithm calculations, use [SciPy](https://www.scipy.org/) and/or [NumPy](https://www.numpy.org/). 
* To perform operations on trees, use [TreeSwift](https://github.com/niemasd/TreeSwift).


## Homework Deliverables
* **[generate_seqs.py](generate_seqs.py):** Your Python 3 code for simulating sequence evolution under the GTR model
    * Usage: `python3 generate_seqs.py -t <tree> -p <gtr_params> -r <root_seq> -o <out_file>`
    * We have provided starter code in the file, but you must fill out the `evolve(tree, root_seq, gtr_probs, gtr_rates)` and `random_seq(k, gtr_probs)` functions
        * **Note:** Feel free to add import statements to NumPy and SciPy within these functions as you need!
* **[writeup.pdf](writeup.pdf):** A write-up that includes the following pieces of information:
    1. A description of how you normalize the GTR transition rates
    2. A description of the time complexity of your code

## Grade Breakdown (100 Points)
* **[writeup.pdf](writeup.pdf): 30 Points**
    * *Valid normalization of GTR transition rates:* 20 points
    * *Correct time complexity:* 10 points

* **[generate_seqs.py](generate_seqs.py): 70 Points**
    * *Code Execution:* 30 points
    * *Short Sequences (k = 200):* 20 points
        * 10 replicates, each worth 2 points
    * *Long Sequences (k = 2,000):* 20 points
        * 10 replicates, each worth 2 points

## References
1. Tavaré S. Some probabilistic and statistical problems in the analysis of DNA sequences. *Lectures on Mathematics in the Life Sciences*, 17:57–86, 1986.
