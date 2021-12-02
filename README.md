# ![seq-toolkit.logo](/docs/images/seq.toolkit.logo.svg)

Basic functions for genomic sequence manipulation

[![PyPI pyversions](https://img.shields.io/pypi/pyversions/seq-toolkit.svg)](https://pypi.python.org/pypi/seq-toolkit/)
[![PyPI version](https://badge.fury.io/py/seq-toolkit.svg)](https://badge.fury.io/py/seq-toolkit)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

### Example use-cases

#### Create a sequence
```python
import seq_toolkit

sequence = seq_toolkit.Seq()
sequence.simulate(n_bases=250, return_seq=True)
```
>'GTTCAACTTTAACCAGCGGTTATGCTCTTCCATGTGAGTCGTAGTCGGGTTCGCACGAAAGAATTATTTATTGCTGCGATACGCAATGCATCTGGTTGTGGAGTTTCACCAGGCAGGCAATTAGTCCTATGCGGAACCTGCTGCTATAAAACGCATAAATTAACTGGCACACCAGGGAGGTAAGGGATGAGAGGCCTACAAGATTCCCATGTGCATATGGAGGGCGTACTGGATTCACGCGTTCGGAGGC'

#### Simulate a gene
```python
import seq_toolkit

gene = seq_toolkit.Gene()
gene.create(n_exons=8, return_gene=False)
gene.plot()
```
![image](https://user-images.githubusercontent.com/47393421/144492953-81b016b7-710e-414f-9e6c-b3576e2fc33c.png)

* For more on ene simulation, see the following: [example notebook](/docs/notebooks/01.GeneSimulation.ipynb)

#### Manipulate a sequence
```python
import seq_toolkit

seq = 'CGGGTTCGCACGAAAGAATTATTTATTGCTGCGATACGCAATGCATCTGGTTGTGGAGTTTCACCAGGCAGGCAATTAGTCCTATGCGGAACCTGCTGCTATAAAACGCAT'
sequence = seq_toolkit.SequenceManipulator(seq)

sequence.reverse_complement()
```
>'ATGCGTTTTATAGCAGCAGGTTCCGCATAGGACTAATTGCCTGCCTGGTGAAACTCCACAACCAGATGCATTGCGTATCGCAGCAATAAATAATTCTTTCGTGCGAACCCG'

### Installation

```python
pip install seq-toolkit
```

To install the development version:

```BASH
git clone https://github.com/mvinyard/seq-toolkit.git

cd ./seq-toolkit/; pip install -e .
```

### Notes
For questions, please send an email to [mvinyard@broadinstitute.org](mailto:mvinyard@broadinstitute.org) or open an [issue](https://github.com/mvinyard/seq-toolkit/issues). 
