# msa-ga
Multiple Sequence Alignment - Genetic Algorithm

This implementation requires **Python 3**. To run, do in your command line:
```
git clone https://github.com/filipefalcaos/msa-ga.git
cd msa-ga
python3 src/main.py -i "data/data_1.txt" -c 100 -g 400 -min 200 -mut 0.03 -p 1.5
```
The above example runs the Genetic Algorithm for the sequences in "data/data_1.txt", using 100 chromosomes, at most 400 generations, a minimum number of 200 generations before stopping, a mutation rate of 0.03, and a penalty of 1.5 for columns that contain only gaps.

**Disclaimer**: This code is not intended to be performatic!
