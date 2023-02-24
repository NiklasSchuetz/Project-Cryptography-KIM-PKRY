#   KIM-PKRY Niklas Schütz 3761908
#   Task 3: Choose k out of n combinations 

from time import perf_counter

def n_ueber_k(n, k):
    ''' function to calculate all combinations of n over k.
    '''
    cur_combination = list(range(k))
    # to not get combinations multiple times we use a set in which we add the combination as a tuple
    solution_set = {tuple(cur_combination)}
    while 1:
        ok = False
        j = 0
        for j in range(k):
            if cur_combination[k-1-j] < (n-1-j):
                cur_combination[k-1-j] += 1
                solution_set.add(tuple(cur_combination))
                ok = True
                break

        if not ok:
            return solution_set

        j -= 1

        while j >= 0:
            cur_combination[k-1-j] = cur_combination[k-2-j] + 1
            solution_set.add(tuple(cur_combination))
            j -= 1

    return solution_set


if __name__ == "__main__":

    n = 10
    k = 3

    start = perf_counter()
    combinations = n_ueber_k(n, k)
    end = perf_counter()

    print(f'{len(combinations)} combinations in: {format(end-start, ".5f")} seconds')
    
    for combination in sorted(combinations):
        print(combination)
    