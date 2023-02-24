# KIM-PKRY HTW SAAR 2021/2022
# Niklas Schütz

from math import log2, ceil, sqrt
from random import sample
from time import perf_counter
import itertools
from collections import defaultdict
from typing import List


# HILFSFUNKTIONEN
def k_aus_n(n, k):
    """function to calculate all combinations of n over k.
    based on pseudocode given in the lecture
    """
    cur_combination = list(range(k))
    # to not get combinations multiple times we use a set in which we add the combination as a tuple
    solution_set = {tuple(cur_combination)}
    while 1:
        ok = False
        j = 0
        for j in range(k):
            if cur_combination[k - 1 - j] < (n - 1 - j):
                cur_combination[k - 1 - j] += 1
                solution_set.add(tuple(cur_combination))
                ok = True
                break

        if not ok:
            return solution_set

        j -= 1

        while j >= 0:
            cur_combination[k - 1 - j] = cur_combination[k - 2 - j] + 1
            solution_set.add(tuple(cur_combination))
            j -= 1

    return solution_set


def create_help_arrays(g, p, max: int):
    set_arr = {}
    unset_arr = {}
    for i in range(0, max):
        set_arr[i] = g ^ (2 ^ i) % p
        unset_arr[i] = g ^ (p - 1 - 2 ^ i) % p

    return set_arr, unset_arr


def get_bit_position(nr: str):
    nr = nr[2:]
    length = len(nr)

    return length - nr.find("1") - 1


def get_bit_positions(nr: str):
    nr = nr[2:]
    length = len(nr)

    return [length - pos - 1 for pos, bit in enumerate(nr) if bit == "1"]


def create_number_of_length_n_with_x_bits(n, x):
    number = 0
    for bit in sample(list(range(n)), x):
        number |= 1 << bit
    return number


def log(str: str, file=None):

    if file != None:
        try:
            with open(f"./logs/{file}.txt", "a") as f:
                f.write(f"{str}\n")
        except FileNotFoundError:
            with open(f"./logs/{file}.txt", "w") as f:
                f.write(f"{str}\n")

    print(str)


def timing_logging_decorator(fnc):
    """decorator to time the function and log the result to stdout"""

    global x_bit_length

    if "x_bit_length" not in globals():
        x_bit_length = "no bitlength given as global value"

    def inner(*args, **kwargs):
        start = perf_counter()
        result = fnc(*args, **kwargs)
        end = perf_counter()

        time = float(format(end - start, ".5f"))

        # log(
        #     # f'{fnc.__name__}: {"correct" if result == x else "no result found"}: {result} in {format(end-start, ".5f")} seconds'
        #     f"{fnc.__name__}: {x_bit_length} bit in {time} seconds {"with"+  kwargs["max_diff"] + "diff" if fnc.__name__ == "dlog_split_combination_best" else ""}",
        #     file=fnc.__name__,
        # )

        if fnc.__name__ == "dlog_split_combination_best":
            print(
                f'{fnc.__name__}: {"correct" if result == x else "no result found"}: {result} in {format(end-start, ".5f")} seconds with max_diff: {kwargs["max_diff"]}'
            )
        else:
            print(
                f'{fnc.__name__}: {"correct" if result == x else "no result found"}: {result} in {format(end-start, ".5f")} seconds'
            )

        return result, time

    return inner


# This function generates all n bit Gray
# codes and prints the generated codes
# This code is contributed
# by Mohit kumar 29
def generateGrayarr(n):

    # base case
    if n <= 0:
        return

    # 'arr' will store all generated codes
    arr = list()

    # start with one-bit pattern
    arr.append("0")
    arr.append("1")

    # Every iteration of this loop generates
    # 2*i codes from previously generated i codes.
    i = 2
    j = 0
    while True:

        if i >= 1 << n:
            break

        # Enter the prviously generated codes
        # again in arr[] in reverse order.
        # Nor arr[] has double number of codes.
        for j in range(i - 1, -1, -1):
            arr.append(arr[j])

        # append 0 to the first half
        for j in range(i):
            arr[j] = "0" + arr[j]

        # append 1 to the second half
        for j in range(i, 2 * i):
            arr[j] = "1" + arr[j]
        i = i << 1

    return arr


# DLOG ALGORITHMEN
@timing_logging_decorator
def babystep_giantstep(g: int, h: int, p: int) -> int:
    """shanks babystep giantstep algorithm to calculate x for a given g^x = h mod p
    based on the pseudocode given in the lecture
    """
    m = ceil(sqrt(p))

    babysteps = {}

    # bs: for i in range 0..m calculate (h*g^p-i-1) mod p and store in a dict (result_for_i : i)
    for i in range(m):
        y = (h * pow(g, p - i - 1, p)) % p
        babysteps[y] = i

    # gs: for j in range 0..m calculate (g^m*j) mod p and check if result is in babystep dict. if hit: x = ((h*g^p-i-1) mod p) + m*j
    for j in range(m):
        y2 = pow(g, m * j, p)

        if y2 in babysteps:
            return babysteps[y2] + m * j


@timing_logging_decorator
def dlog_combinations(t, g, h, p, n=None):
    """calculating g^x mod p for all possible x assuming we now how many bits are set in x"""
    hp = h % p

    if n == None:
        n = int(log2(p - 1))

    for s in k_aus_n(n, t):  # iterieren über alle kombinationen
        possible_x = 0  #  kombination in zahl umwandeln
        for bit_pos in s:
            possible_x = possible_x | 1 << bit_pos

        if pow(g, possible_x, p) == hp:  # test ob ergebnis
            return possible_x


@timing_logging_decorator
def dlog_gray_codes(t, g, h, p, n=None):
    """iterating over all graycodes and calculating the  nth result by multiplying the result of the n-1th result with a value from a help array
    two help arrays: one for setting a bit, one for unsetting
    """
    hp = h % p

    if n == None:
        n = int(log2(p - 1))

    unset_arr: List[int] = [
        pow(g, (p - 1 - (2 ** (i - 1))), p) for i in range(1, n + 1)
    ]
    set_arr: List[int] = [pow(g, 2 ** (i - 1), p) for i in range(1, n + 1)]

    graycodes = generateGrayarr(n)

    old_code = int(graycodes.pop(0), 2)
    old_code_res = pow(g, old_code, p)

    if old_code_res == hp:
        return int(old_code)

    for code in graycodes:
        code = int(code, 2)
        diff = code ^ old_code  # get difference
        pos = get_bit_position(
            bin(diff)
        )  # get position of bit which is not the same in current and old graycode

        if diff & code:
            code_res = (old_code_res * set_arr[pos]) % p
        else:
            code_res = (old_code_res * unset_arr[pos]) % p

        if code_res == hp:
            return code

        old_code = code
        old_code_res = code_res


@timing_logging_decorator
def dlog_shanks_mit_splitting(t, g, h, p, n=None):
    """ "normal" shanks but with optimizing by splitting the combinations"""
    if n == None:
        n = int(log2(p - 1))

    # arr_0 = itertools.combinations(range(n), t // 2)
    # arr_1 = itertools.combinations(range(n), t + 1 // 2)

    arr_0 = itertools.combinations(range((n // 2) + 1), t // 2)
    arr_1 = itertools.combinations(range((n // 2) - 1, n), (t + 1) // 2)

    babysteps = {}

    g_inverted = pow(g, -1, p)  # its faster to calculate then g^-1^x1  g^-x1

    for s in arr_0:
        possible_x = 0
        for bit_pos in s:
            possible_x = possible_x | 1 << bit_pos

        babysteps[pow(g, possible_x, p)] = possible_x

    for s in arr_1:
        possible_x = 0
        for bit_pos in s:
            possible_x = possible_x | 1 << bit_pos

        if ((h * pow(g_inverted, possible_x, p)) % p) in babysteps:
            return babysteps[((h * pow(g_inverted, possible_x, p)) % p)] + possible_x


@timing_logging_decorator
def dlog_combinations_with_set_unset_array(t, g, h, p, n=None):
    """IDEA:
    * use n over t/2 splitting
    * like with graycodes we calculate the result by multiplying the old result with one of the help arrays
    """

    if n == None:
        n = int(log2(p - 1))

    arr_0 = itertools.combinations(range(n), t // 2)
    arr_1 = itertools.combinations(range(n), t + 1 // 2)

    # arr_0 = itertools.combinations(range((n // 2) + 1), t // 2)
    # arr_1 = itertools.combinations(range((n // 2) - 1, n), (t + 1) // 2)

    babysteps = {}

    g_inverted = pow(g, -1, p)  # its faster to calculate then g^-1^x1  g^-x1

    set0 = [pow(g, pow(2, i), p) for i in range(0, n)]
    unset0 = [pow(g, (p - 1 - 2 ** (i)), p) for i in range(0, n)]
    set1 = [pow(g_inverted, pow(2, i), p) for i in range(0, n)]
    unset1 = [pow(g_inverted, (p - 1 - 2 ** (i)), p) for i in range(0, n)]

    old_x0, old_x1 = 0, 0
    old_x0_res, old_x1_res = 1, 1

    for combination in arr_0:
        possible_x = 0
        x0_res = old_x0_res

        for i in combination:
            possible_x = possible_x | 1 << i

        bit_positions = get_bit_positions(bin(possible_x ^ old_x0))

        for bit_pos in bit_positions:
            if (1 << bit_pos) & possible_x:
                x0_res *= set0[bit_pos]
            else:
                x0_res *= unset0[bit_pos]
            x0_res %= p

        babysteps[x0_res] = possible_x

        old_x0 = possible_x
        old_x0_res = x0_res

    for s in arr_1:
        possible_x = 0
        x1_res = old_x1_res

        for bit_pos in s:
            possible_x = possible_x | 1 << bit_pos

        bit_positions = get_bit_positions(bin(possible_x ^ old_x1))

        for bit_pos in bit_positions:
            if 1 << bit_pos & possible_x:
                x1_res *= set1[bit_pos]
            else:
                x1_res *= unset1[bit_pos]
            x1_res %= p

        if (h * x1_res) % p in babysteps:
            return babysteps[(h * x1_res) % p] + possible_x

        old_x1 = possible_x
        old_x1_res = x1_res

    return None


@timing_logging_decorator
def dlog_split_combination_diff(t, g, h, p, n=None, max_diff=0):
    """IDEA:
    * use n over t/2 splitting
    * like with graycodes we calculate the result by multiplying the old result with one of the help arrays
    """

    if n == None:
        n = int(log2(p - 1))

    arr_0 = itertools.combinations(range((n // 2) + 1), t // 2)
    arr_1 = itertools.combinations(range((n // 2) - 1, n), (t + 1) // 2)

    # arr_0 = itertools.combinations(range(n), t // 2)
    # arr_1 = itertools.combinations(range(n), t + 1 // 2)

    babysteps = {}

    g_inverted = pow(g, -1, p)  # its faster to calculate then g^-1^x1  g^-x1

    set0 = [pow(g, pow(2, i), p) for i in range(0, n)]
    unset0 = [pow(g, (p - 1 - 2 ** (i)), p) for i in range(0, n)]
    set1 = [pow(g_inverted, pow(2, i), p) for i in range(0, n)]
    unset1 = [pow(g_inverted, (p - 1 - 2 ** (i)), p) for i in range(0, n)]

    old_x0, old_x1 = 0, 0
    old_x0_res, old_x1_res = 1, 1

    old_combination = set()

    for combination in arr_0:
        possible_x = 0
        x0_res = old_x0_res

        for i in combination:
            possible_x = possible_x | 1 << i

        if (
            len(set(combination).symmetric_difference(set(old_combination))) < max_diff
        ):  # check if difference in old and cur comb is smaller than max diff -> use array
            bit_positions = get_bit_positions(bin(possible_x ^ old_x0))

            for (
                bit_pos
            ) in bit_positions:  # iterate over changed bits and set/unset them
                if (1 << bit_pos) & possible_x:
                    x0_res *= set0[bit_pos]
                else:
                    x0_res *= unset0[bit_pos]
                x0_res %= p

        else:
            x0_res = pow(g, possible_x, p)

        babysteps[x0_res] = possible_x

        old_x0 = possible_x
        old_x0_res = x0_res
        old_combination = combination

    old_combination = set()
    for combination in arr_1:

        possible_x = 0
        x1_res = old_x1_res

        for bit_pos in combination:
            possible_x = possible_x | 1 << bit_pos

        if len(set(combination).symmetric_difference(set(old_combination))) < max_diff:
            bit_positions = get_bit_positions(bin(possible_x ^ old_x1))

            for bit_pos in bit_positions:
                if 1 << bit_pos & possible_x:
                    x1_res *= set1[bit_pos]
                else:
                    x1_res *= unset1[bit_pos]
                x1_res %= p

        else:
            x1_res = pow(g_inverted, possible_x, p)

        if (h * x1_res) % p in babysteps:  # check if in babystep
            return babysteps[(h * x1_res) % p] + possible_x

        old_x1 = possible_x
        old_x1_res = x1_res
        old_combination = combination

    return None


if __name__ == "__main__":
    global x_bit_length

    combination_results = defaultdict(list)
    shanks_split_results = defaultdict(list)
    gray_codes_results = defaultdict(list)
    comb_set_unset_results = defaultdict(list)

    p_bit_dict = {
        20: 568979,
        21: 1093739,
        22: 2142299,
        23: 4240043,
        24: 12627779,
        25: 21015803,
        26: 37793219,
        27: 71347643,
        28: 205566419,
        29: 339784427,
        30: 608218643,
        31: 1681961363,
        32: 2755703579,
        33: 4903188827,
        34: 13493121347,
        35: 22083055859,
        36: 56442793907,
        37: 90802537283,
        38: 228241493603,
        39: 365680439459,
        40: 915436254683,
        41: 1465192069163,
        42: 2564703694403,
        43: 6962750206547,
        44: 15758843227979,
        45: 26491917330443,
        46: 44084103375827,
        47: 114452847557219,
        48: 255190335908699,
        49: 536665312624307,
        50: 818140289331083,
        51: 1944040196172827,
        52: 4195840009864523,
        53: 6447639823550603,
        54: 15454839078283427,
        55: 33469237587774683,
        56: 69498034606728779,
        57: 105526831625695619,
        58: 249642019701548819,
        59: 537872395853268563,
        60: 826102772004977387,
        61: 1979024276611826603,
        62: 3131945781218668187,
        63: 7743631799646055499,
        64: 12355317818073449027,
    }

    # set p and g
    p_bit = 48
    p = p_bit_dict[p_bit]
    g = 111

    for x_bit_length in range(20, 21):

        for run in range(1, 2):

            x = create_number_of_length_n_with_x_bits(n=p_bit, x=x_bit_length)
            # 56590574651753 n 48 t 20
            t = x_bit_length
            n = p_bit

            print(f"n: {n} t: {t}")

            # calculate h
            h = pow(g, x, p)

            _, t_com_set_unset = dlog_combinations_with_set_unset_array(
                t=t, g=g, h=h, p=p, n=n
            )
            _, t_shanks_split = dlog_shanks_mit_splitting(t=t, g=g, h=h, p=p, n=n)

            _, t_2 = dlog_split_combination_diff(
                t=t, g=g, h=h, p=p, n=n, max_diff=t // 2
            )
            _, t_4 = dlog_split_combination_diff(
                t=t, g=g, h=h, p=p, n=n, max_diff=t // 4
            )
            _, t_0 = dlog_split_combination_diff(t=t, g=g, h=h, p=p, n=n, max_diff=0)
            _, t_full = dlog_split_combination_diff(
                t=t, g=g, h=h, p=p, n=n, max_diff=10000
            )

        #     # babystep_giantstep(g=g, h=h, p=p)
        #     _, t_comb = dlog_combinations(t=t, g=g, h=h, p=p, n=n)
        #     _, t_gray = dlog_gray_codes(t=t, g=g, h=h, p=p, n=n)
