import time
import argparse
from ga import GA


def main():

    # Create the input parser
    parser = argparse.ArgumentParser(description="Mine pull requests from GitHub's open repositories.")

    parser.add_argument("-input", metavar='-i', type=str, help="Input path")
    parser.add_argument("-chromosomes", metavar='-c', type=int, help="Number of chromosomes")
    parser.add_argument("-gens", metavar='-gen', type=int, help="Number of generations")
    parser.add_argument("-min_gens", metavar='-min', type=int, help="Minimum number of generations before stopping")
    parser.add_argument("-mutation_rate", metavar='-mut', type=float, help="Mutation rate")

    # Parse args
    args = parser.parse_args()

    # Run the GA
    start = time.time()
    genetic_algorithm = GA(args.chromosomes, args.gens, args.min_gens, args.mutation_rate)
    genetic_algorithm.run_ga(args.input)
    end = time.time()
    print("\nRunning time: " + str(end - start) + " seconds")


if __name__ == "__main__":
    main()
