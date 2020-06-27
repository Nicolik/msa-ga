# Class that implements the proposed genetic algorithm
# The initial population, selection, substitution, adaptation function,
# crossover and mutation operators are implemented in this class.

import nwalign3 as nw
import numpy as np
import random
import math
from utils import Utils


class GA:

    def __init__(self, chromosomes, generations, min_generations, mutation_rate):
        self.chromosomes = chromosomes
        self.generations = generations
        self.min_generations = min_generations
        self.mutation_rate = mutation_rate

    # Use nwalign to compute the score of the pairwise alignments
    # using a BLOSUM62 scoring matrix
    @staticmethod
    def evaluation_func(lines):
        for i in range(len(lines)):
            sum_score = 0

            # Compute the pairwise scores
            for j in range(len(lines)):
                if i != j:
                    scr = nw.score_alignment(lines[i], lines[j], gap_open=-1, gap_extend=-0.5, matrix='matrix/BLOSUM62')
                    sum_score += scr

            return sum_score

    # Returns if there are no relevant changes
    # after 100 generations
    def no_change(self, best):
        if len(best) < self.min_generations:
            return False
        else:
            percent = int(0.2 * len(best))
            last = best[-percent:]

            if np.var(last) < 1.05:
                return True

        return False

    # Create the initial population
    def init_pop(self, lines_list):
        pop = []
        for c in range(self.chromosomes):

            # New chromosome
            lines_list_aux = []

            # Use nwalign to compute the pairwise alignments
            # by the Needleman-Wunsch algorithm
            for i in range(len(lines_list)):
                alignments = []

                # Compute the pairwise alignments
                for j in range(len(lines_list)):
                    if i != j:
                        curr_alignment = nw.global_align(lines_list[i], lines_list[j])
                        alignments.append(curr_alignment)

                # Randomly select an alignment
                alignment = random.choice(alignments)
                alignment = alignment[0]
                lines_list_aux.append(alignment)

            # Add the chromosome generated and prints it
            lines_list_aux = Utils.add_gaps(lines_list_aux)
            pop.append({"chromosome": lines_list_aux, "evaluation": 0})
            print("\nChromosome " + str(c + 1) + ":")
            Utils.print_chromosome(lines_list_aux)

        # Initial population
        return pop

    # Select two parents for a crossover operation
    @staticmethod
    def select_parents(pop, evaluations):

        # Normalize the fitness
        pop_sum = sum(evaluations)
        evaluations_aux = []
        for ev in evaluations:
            evaluations_aux.append(ev / pop_sum)

        # Build the mating pool
        pool = []
        for i in range(len(pop)):
            prob = math.ceil(evaluations_aux[i] * 100)

            for j in range(prob):
                pool.append(pop[i])

        # Randomly select two parents
        p1 = random.choice(pool)
        p2 = random.choice(pool)
        return [p1, p2]

    # Apply a crossover operation on p1 and p2
    @staticmethod
    def apply_crossover(pop, p1, p2):

        n = len(pop[0]["chromosome"])

        # Apply horizontal crossover
        rand_h = random.randint(1, n - 1)
        # print("\nSlice horizontally by: " + str(rand_h))
        child = p1["chromosome"][:rand_h] + p2["chromosome"][rand_h:]

        # Return the child
        return {"chromosome": child, "evaluation": 0}

    # Apply a mutation on child
    def apply_mutation(self, pop, child):

        n = len(pop[0])

        # Select mutation method and apply
        if round(random.uniform(0, 1), 2) < self.mutation_rate:
            rand = round(random.uniform(0, 1), 2)

            # Apply gaps removal
            if rand < 0.5:
                cell_i = random.randint(1, n - 1)
                cell_j = random.randint(1, len(child[cell_i]) - 1)

                while child[cell_i][cell_j] != "-":
                    cell_i = random.randint(0, n - 1)
                    cell_j = random.randint(0, len(child[cell_i]) - 1)

                # Get extending gap
                gaps = Utils.get_interval_gaps(child, cell_i, cell_j)
                start, end = gaps[0], gaps[len(gaps) - 1]

                # Remove gaps
                child[cell_i] = child[cell_i][:start] + child[cell_i][end + 1:]

                # if start == end:
                #     print("\nMutate: removed gap from (" + str(cell_i + 1) + ", " + str(start + 1) + ")")
                # else:
                #     print("\nMutate: removed gaps from (" + str(cell_i + 1) + ", " + str(start + 1) + ") until (" +
                #           str(cell_i + 1) + ", " + str(end + 1) + ")")

            # Apply k addition
            else:
                cell_i = random.randint(1, n - 1)
                cell_j = random.randint(1, len(child[cell_i]) - 1)
                k = random.randint(1, math.ceil(0.1 * Utils.calc_m(child)))
                to_add = ""

                for i in range(k):
                    to_add += "-"

                child[cell_i] = child[cell_i][:cell_j] + to_add + child[cell_i][cell_j:]
                # print("\nMutate: added " + str(k) + " gaps in (" + str(cell_i + 1) + ", " + str(cell_j + 1) + ")")

        # Return the child
        return child

    def run_ga(self, input_path):

        # Read input file string
        lines_list = Utils.prepare_input(input_path)

        # Prints the original chromosome
        print("Input matrix:")
        Utils.print_chromosome(lines_list)

        # Create the initial population
        pop = self.init_pop(lines_list)

        # Repeat for all generations or until a good solution
        # appears
        best_val = None
        best_chromosome = None
        best_chromosomes = []
        count = 0
        new_pop = []
        print()

        while count < self.generations:

            # Evaluate all chromosomes
            evaluations = []
            for i in range(len(pop)):
                pop[i]["evaluation"] = self.evaluation_func(pop[i]["chromosome"])
                evaluations.append(pop[i]["evaluation"])

            # Repeat for all chromosomes
            for w in range(self.chromosomes):

                # Select two parents
                p1, p2 = self.select_parents(pop, evaluations)
                # print("\nParent 1:")
                # Utils.print_chromosome(p1["chromosome"])
                # print("\nParent 2:")
                # Utils.print_chromosome(p2["chromosome"])

                # Get crossover probability
                rand = round(random.uniform(0, 1), 2)

                # Apply a crossover operation on p1 and p2
                if rand < 0.5:
                    child = self.apply_crossover(pop, p1, p2)

                    # Apply a mutation on child
                    child["chromosome"] = self.apply_mutation(pop, child["chromosome"])

                    # Add the child to the new population
                    child["chromosome"] = Utils.add_gaps(child["chromosome"])
                    child["chromosome"] = Utils.remove_useless_gaps(child["chromosome"])
                    new_pop.append(child)
                    # print("\nChild:")
                    # Utils.print_chromosome(child["chromosome"])

            # Get the best chromosome
            best_val = 0
            best_chromosome = None
            for i in range(len(new_pop)):
                curr_val = self.evaluation_func(new_pop[i]["chromosome"])
                new_pop[i]["evaluation"] = curr_val

                if curr_val >= best_val:
                    best_val = curr_val
                    best_chromosome = new_pop[i]["chromosome"]

            # Print stats
            print("Generation " + str(count + 1) + ": " + str(best_val))

            # Add best to the list
            best_chromosomes.append(best_val)

            # Break the execution if there are no
            # relevant changes
            if self.no_change(best_chromosomes):
                print("\nAbort: no variation!", end="")
                break

            # Update population
            # Add to the population the new generated chromosomes and
            # remove the same number of the worst chromosomes from
            # the original population
            pop = sorted(pop, key=lambda k: k["evaluation"], reverse=True)
            new_pop = sorted(new_pop, key=lambda k: k["evaluation"], reverse=True)

            for i in range(len(new_pop)):
                pop.pop()
                pop.insert(0, new_pop[i])

            new_pop = []
            count = count + 1

        # Best solution
        print("\nBest solution:")
        Utils.print_chromosome(best_chromosome)
        return best_val
