# Class that handles utility methods
#
# More in detail, it handles mostly the addition, counting and
# removal of gaps and spaces from a matrix. It also reads and prepares
# the input matrix for the GA.

import os
import sys
from contextlib import contextmanager


class Utils:

    # Prints a chromosome's sequences
    @staticmethod
    def print_chromosome(chromosome):
        for seq in chromosome:
            print(seq)

    # Calc the maximum number of columns m
    @staticmethod
    def calc_m(lines):
        m_aux = 0
        lengths = []

        for line in lines:
            length = len(line)
            lengths.append(length)

            if length >= m_aux:
                m_aux = length

        return m_aux

    # Add the required gaps to produce an
    # initial alignment
    @staticmethod
    def add_gaps(lines):

        # Complete with gaps
        for i in range(len(lines)):
            diff = Utils.calc_m(lines) - len(lines[i])
            for j in range(diff):
                lines[i] += "-"

        return lines

    # Remove useless gaps from an alignment
    @staticmethod
    def remove_useless_gaps(lines):
        to_rm = []
        only_gaps = True
        m = Utils.calc_m(lines)

        # Find columns with gaps only
        for i in range(m):
            for j in range(len(lines)):
                if (lines[j][i] != " ") and (lines[j][i] != "-"):
                    only_gaps = False

            if only_gaps:
                to_rm.append(i)
            else:
                only_gaps = True

        # Remove gap-only columns
        to_rm.reverse()
        for i in range(len(lines)):
            line = []
            for j in range(m):
                line.append(lines[i][j])

            for j in to_rm:
                line.pop(j)

            lines[i] = ''.join(line)

        # Return the new matrix
        return lines

    # Get an interval of gaps by a given gap
    @staticmethod
    def get_interval_gaps(lines, cell_i, cell_j):
        aux_cell_j = cell_j
        gaps = []
        symbol_found = False

        # Find gaps after cell_j (inclusive)
        while not symbol_found:
            if cell_j > len(lines[cell_i]) - 1:
                break
            elif lines[cell_i][cell_j] == "-":
                gaps.append(cell_j)
                cell_j += 1
            else:
                symbol_found = True

        # Find gaps before cell_j (exclusive)
        aux_cell_j -= 1
        while not symbol_found:
            if aux_cell_j < 0:
                break
            elif lines[cell_i][aux_cell_j] == "-":
                gaps.insert(0, aux_cell_j)
                aux_cell_j -= 1
            else:
                symbol_found = True

        return gaps

    @staticmethod
    def prepare_input(input_path):

        # Read input file string
        with open(input_path) as f:
            input_str = f.read()

        # Get lines
        input_str = input_str.replace(" ", "")
        lines_list = input_str.split("\n")

        # Add the gaps
        lines_list = Utils.add_gaps(lines_list)
        return lines_list

    @staticmethod
    @contextmanager
    def suppress_stdout():
        with open(os.devnull, "w") as devnull:
            old_stdout = sys.stdout
            sys.stdout = devnull
            try:
                yield
            finally:
                sys.stdout = old_stdout
