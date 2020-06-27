from utils import Utils
from ga import GA


# ClustalW output
clustal_1 = Utils.prepare_input("data/sequences/clustal_2.out")
genetic_algorithm = GA(None, None, None, None)
print(genetic_algorithm.evaluation_func(clustal_1))
