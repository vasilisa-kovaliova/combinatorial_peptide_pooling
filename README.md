# Combinatorial peptide pooling

## Intro
A script for the distribution of peptides between N_POOLS number of pools where each peptide is added to ITER number of pools.

### Task
Let's imagine that we want to find an epitope within some protein that is recognized by the TCR of interest. To do so, we need to cut this protein into multiple peptides and test a T-cell bearing a TCR of interest to find the epitope recognized by this TCR. To be more efficient, we can distribute these peptides across some appropriate amount of pools, and this way, we will reduce the number of tests.

### Example 1: Non-overlapping peptides
There are 6 peptides. We want to distribute them across 4 pools (N_POOLS = 4), and one peptide is added to 2 pools (ITER = 2)
Here is a possible scheme:

| Peptide           | Pool 1 | Pool 2 | Pool 3 | Pool 4 |
| ----------------- | ------ | ------ | ------ | ------ |
| MFVFLVLLPLVSSQCVN | X      | X      |        |        |
| VYYPDKVFRSSVLHSTQ | X      |        | X      |        |
| NPVLPFNDGVYFASTEK | X      |        |        | X      |
| FGTTLDSKTQSLLIVNN |        | X      | X      |        |
| LGVYYHKNNKSWMESEF |        | X      |        | X      |
| VSQPFLMDLEGKQGNFK |        |        | X      | X      |

This way, after the co-culture of these pools with T-cells bearing a TCR  of interest, we can understand which peptide led to the activation of this T-cell.

| Pools activated    | Peptide           | 
| ------------------ | ----------------- |
| 1 and 2            | MFVFLVLLPLVSSQCVN |
| 1 and 3            | VYYPDKVFRSSVLHSTQ |
| 1 and 4            | NPVLPFNDGVYFASTEK |
| 2 and 3            | FGTTLDSKTQSLLIVNN |
| 2 and 4            | LGVYYHKNNKSWMESEF |
| 3 and 4            | VSQPFLMDLEGKQGNFK |

### Example 2: Overlapping peptides
But what if the cognate epitope of tested TCR is located on the junction of two peptides? To rule out this option, we need to use overlapping peptides with appropriate overlap length. Then we need to distribute them again across N_POOLS number of pools.
For example like this:

| Peptide           | Pool 1 | Pool 2 | Pool 3 | Pool 4 |
| ----------------- | ------ | ------ | ------ | ------ |
| MFVFLVLLPLVSSQCVN | X      | X      |        |        |
| VLLPLVSSQCVNLTTRT | X      |        | X      |        |
| VSSQCVNLTTRTQLPPA | X      |        |        | X      |
| VNLTTRTQLPPAYTNSF |        | X      | X      |        |
| RTQLPPAYTNSFTRGVY |        | X      |        | X      |
| PAYTNSFTRGVYYPDKV |        |        | X      | X      |

If a recognized epitope is located in the overlapping part of two peptides, it may lead to confusing results:
| Right Epitope  | Pools activated | Peptide                               |
| ---------------| --------------- | ------------------------------------- |
| VLLPLVSS       | 1, 2, 3         | MFVFLVLLPLVSSQCVN or VLLPLVSSQCVNLTTRT|
| VSSQCVNL       | 1, 3, 4         | VLLPLVSSQCVNLTTRT or VSSQCVNLTTRTQLPPA|
| VNLTTRTQ       | 1, 2, 3, 4      | VSSQCVNLTTRTQLPPA or VNLTTRTQLPPAYTNSF|
| RTQLPPAY       | 2, 3, 4         | VNLTTRTQLPPAYTNSF or RTQLPPAYTNSFTRGVY or PAYTNSFTRGVYYPDKV|

In the first 3 rows, we still can understand what sequence led to TCR activation, but in the case of the last two rows, we can not distinguish between VNLTTRTQLPPAYTNSF, RTQLPPAYTNSFTRGVY, and PAYTNSFTRGVYYPDKV peptides.

To avoid such situations, overlapping peptides should be distributed across pools carefully.
Provided script does so.

## Algorithm
To avoid situations like in example 2 (we call them collisions), this script assigns peptides to pools in order. To make number of peptides in each pool more or less balanced, random factor is involved.

After the first round of peptide assignment, check for collisions happens. In case of collisions, peptides are re-assigned until the number of collisions reaches the specified number. If, after multiple reassignments, collisions are present, the script stops.

### Installation
```
git clone https://github.com/vasilisa-kovaliova/combinatorial_peptide_pooling.git
```
### Requirements
Pandas and numpy modules.

You can install them by typing in command line:
```
pip install -r requirements.txt
```

### Usage
python Peptide_pooling.py -n_pools N_POOLS -iters ITERS -peptides PEPTIDES -unresolved MAX_COLLISIONS -pools OUTPUT_FILE_1 -simulation OUTPUT_FILE_2

Parameters:

--  n_pools -- the number of pools

-- iters -- to how many pools one peptide is added

-- peptides -- path to the file with ordered peptides, columns must be separated by tabs, peptides should be in the column named "SequenceAsEntered"

-- unresolved -- maximum number of peptides that can be undistinguished; this number is helpful if you know that some of your peptides overlap by more than 12 amino acids

-- pools -- path to the file where pools will be written

-- simulation -- path to the file where simulation results will be written, with results for each possible epitope of length 8

N_POOLS and ITERS should be chosen carefully to provide enough pool combinations and not have a too high number of peptides in one pool.

### Possible options

| Number of peptides | N_POOLS | ITERS | number of combinations     | approximate number of peptides in one pool |
| ------------------ | ------- | ----- | -------------------------- | ------------------------------------------ |
| 250                | 12      | 4     | combinations(12, 4) = 495  | (250*4)/12 = 83                            |
| 400                | 16      | 3     | combinations(16, 3) = 560  | (400*3)/16 = 75                            |
| 800                | 15      | 4     | combinations(15, 4) = 1365 | (800*4)/15 = 213                           |
| 1500               | 18      | 4     | combinations(18, 4) = 3060 | (1500*4)/18 = 333                          |
| 2000               | 16      | 6     | combinations(16, 6) = 8008 | (2000*6)/16 = 750                          |


### Example
```
python Peptide_pooling.py -n_pools 12 -iters 4 -peptides peptides_example.tsv -unresolved 0 -pools pools.txt -simulation simulation.tsv
```
