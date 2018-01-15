import sys
from kemeny import KemenyRanking


if __name__ == '__main__':
    fp = sys.argv[1]

    kemeny = KemenyRanking(fp)
    kemeny.build_Q()
    kemeny.solve_ilp()
    kemeny.postprocess()
    kemeny.print_sol()
