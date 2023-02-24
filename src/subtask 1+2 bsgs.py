#   KIM-PKRY Niklas Schütz 3761908
#   Task 1+2: Babystep Giant Step 

from math import ceil, sqrt
from time import perf_counter

def babystep_giantstep(g: int,h: int,p: int) -> int:
    '''shanks babystep giantstep algorithm to calculate x for a given g^x = h mod p 
       based on the pseudocode given in the lecture
    '''
    m = ceil(sqrt(p))

    babysteps = {}

    # bs: for i in range 0..m calculate (h*g^p-i-1) mod p and store in a dict (result_for_i : i)
    for i in range(m):
        y = (h*pow(g, p-i-1, p)) % p 
        babysteps[y] = i

    # gs: for j in range 0..m calculate (g^m*j) mod p and check if result is in babystep dict. if hit: x = ((h*g^p-i-1) mod p) + m*j
    for j in range(m):
        y2 = pow(g, m*j, p)

        if y2 in babysteps:
            return babysteps[y2]+m*j


if __name__ == "__main__":
    timing_dict = {
        20 : 568979,
        21 : 1093739,
        22 : 2142299,
        23 : 4240043,
        24 : 12627779,
        25 : 21015803,
        26 : 37793219,
        27 : 71347643,
        28 : 205566419,
        29 : 339784427,
        30 : 608218643,
        31 : 1681961363,
        32 : 2755703579,
        33 : 4903188827,
        34 : 13493121347,
        35 : 22083055859,
        36 : 56442793907,
        37 : 90802537283,
        38 : 228241493603,
        39 : 365680439459,
        40 : 915436254683,
        41 : 1465192069163,
        42 : 2564703694403,
        43 : 6962750206547,
        44 : 15758843227979,
        45 : 26491917330443,
        46 : 44084103375827,
        47 : 114452847557219,
        48 : 255190335908699,
        49 : 536665312624307,
        50 : 818140289331083,
        51 : 1944040196172827,
        52 : 4195840009864523,
        53 : 6447639823550603,
        54 : 15454839078283427,
        55 : 33469237587774683,
        56 : 69498034606728779,
        57 : 105526831625695619,
        58 : 249642019701548819,
        59 : 537872395853268563,
        60 : 826102772004977387,
        61 : 1979024276611826603,
        62 : 3131945781218668187,
        63 : 7743631799646055499,
        64 : 12355317818073449027
    }

    # set g and h to be constant over all iterations
    g = 187
    x = 420

    # iterate over all provided prime numbers and calculate dlog187(h) mod p
    # Write to file if x is correct and how long the script took for every calculation

    for pos, (k, v) in enumerate(timing_dict.items()):
        # set p to the prime number given in the dict and calculate an h from g,x,p so we can easily check if our algorithm calculated the correct x
        p = v
        h = pow(g,x,p)

        start = perf_counter()
        calculated_x = babystep_giantstep(g=g,h=h,p=p)
        end = perf_counter()
        with open("timingpypy.txt","a") as file:
            file.write(f'{k} bit in:\t{format(end-start, ".5f")} seconds\n')
        print(f'{"correct" if x == calculated_x else "wrong"}: {k} bit in {format(end-start, ".5f")} seconds')

