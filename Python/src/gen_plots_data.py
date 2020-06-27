import glob
import numpy as np
from ga import GA
from utils import Utils


def mutation_data():

    # Plot parameters
    data_source = "data/sequences/"
    gens = 450
    chromosomes = 250
    mutations = [0.01, 0.015, 0.02, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05, 0.055, 0.06, 0.065, 0.07, 0.075, 0.08, 0.085,
                 0.09, 0.095, 0.1]

    # Create output file
    with open("data/output_data2.csv", 'w') as output_file:
        output_file.write("mutation, result\n")

    # Get files in the data source
    files = glob.glob(data_source + "*.txt")

    # Generate data
    for mutation in mutations:
        results = []
        print("Running for mutation rate = " + str(mutation))

        for file in files:
            with Utils.suppress_stdout():
                genetic_algorithm = GA(chromosomes, gens, gens, mutation)
                result = genetic_algorithm.run_ga(file)
                results.append(result)

        # Calc median
        median = np.median(results)

        # Write data to file
        with open("data/output_data2.csv", 'a') as output_file:
            output_file.write(str(mutation) + ", " + str(median) + "\n")


def generation_data():

    # Plot parameters
    data_source = "data/sequences/"
    gens = [50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 350, 375, 400, 425, 450, 475, 500, 525]
    chromosomes = [10, 15, 20, 25, 30, 35, 40, 45, 50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300, 325]

    # Create output file
    with open("data/output_data.csv", 'w') as output_file:
        output_file.write("gens, chromosomes, result\n")

    # Get files in the data source
    files = glob.glob(data_source + "*.txt")

    # Generate data
    for i in range(len(gens)):
        results = []
        print("Running for gen = " + str(gens[i]) + " and chromosomes = " + str(chromosomes[i]))

        for file in files:
            with Utils.suppress_stdout():
                genetic_algorithm = GA(chromosomes[i], gens[i], gens[i], 0.05)
                result = genetic_algorithm.run_ga(file)
                results.append(result)

        # Calc median
        median = np.median(results)

        # Write data to file
        with open("data/output_data.csv", 'a') as output_file:
            output_file.write(str(gens[i]) + ", " + str(chromosomes[i]) + ", " + str(median) + "\n")


if __name__ == "__main__":
    # generation_data()
    # print()
    mutation_data()
