#!/usr/bin/env python3

import os
import time
import argparse
import tempfile
import numpy as np
from mpi4py import MPI
from scipy.special import comb

comm = MPI.COMM_WORLD
comm_rank = comm.Get_rank()
comm_size = comm.Get_size()


class information:
    def version(self):
        version = '3.0.1'
        print(version)
    def auther(self):
        auther = 'Ke Zhang'
        print(auther)
    def email(self):
        email = 'zk2973_p@163.com'
        print(email)
    def license(self):
        license = 'GPL v3'
        print(license)
    def description(self):
        description = 'The acceleration ratio of this MPI script is close to one, and it fully supports pypy3. It is suitable for calculating the exhaustive number of combinations on HPC. For example: 400 threads processing comb (6,1000) takes about two to three days, of course, it depends on what you need to calculate.'
        print(description)
    def algorithm(self):
        algorithm = '''My colleague wrote an exhaustive script, but it took too long to run (using python's itertools package). After I took over, I did a thorough efficiency optimization. The specific steps are as follows:
        1.Obtain the combination number information, estimate and divide the combination number according to the number of CPU threads.
        2. Create and distribute message objects needed by all processes.
        3. The process receives the message, initializes the exhaustion process, and then starts the exhaustion.'''
        print(algorithm)

class message:
    def __init__(self, begin, end, m, n):
        self.begin = begin
        self.end = end
        self.m = m
        self.n = n
        self.final = self.__create_final_comb(m, n)
    def __create_final_comb(self, m, n):
        final_list = []
        final_index = 1
        for i in range(0,n):
            final_list.append(m - n + final_index)
            final_index += 1
        return final_list

def update_begin_comb(a, b):
    for i in range(len(a)-1,-1,-1):
        if a[i] < b[i]:
            if i < len(a)-1:
                a[i] += 1
                for h in range(i+1,len(a)):
                    a[h] = a[i] + h - i
                return a
            else:
                a[i] += 1
                return a
        else:
            continue
    return None

def _argparse():
    parser = argparse.ArgumentParser(description='Exhaustive combination MPI program')
    parser.add_argument('--m', type=int, required=True, help='The length of the combination')
    parser.add_argument('--n', type=int, required=True, help='Number of elements')
    return parser.parse_args()

def split_comb(m, n):
    begin_comb = [1,2]
    end_comb = [m-n+1,m-n+2]
    other_digit = n - 2

    total_index = 0
    result_hash = {}
    while (not begin_comb == end_comb):
        total_index += comb(m-begin_comb[1],other_digit)
        result_hash[total_index] = begin_comb[:]
        #update begin_comb
        begin_comb = update_begin_comb(begin_comb, end_comb)
    #print(" ".join(str(i) for i in begin_comb))
    total_index += comb(m-begin_comb[1],other_digit)
    result_hash[total_index] = begin_comb[:]
    return result_hash

def get_message(drug_number, comb_digit, comm_size):
    split_hash = split_comb(drug_number,comb_digit)

    comb_number = comb(drug_number,comb_digit)
    split_array = np.linspace(0,comb_number,comm_size+1)
    split_array = list(map(int, np.delete(split_array,0)))
    message_array = []
    times_number = 0
    begin_site = [1,2]
    for key,value in split_hash.items():
        #print("key: %d,value: %s" % (key,value))
        if key >= split_array[times_number]:
            times_number += 1
            every_message = message(begin_site[:],value[:],drug_number,comb_digit)
            begin_site = value[:]
            message_array.append(every_message)
        else:
            continue
    return message_array

def compound():
    i_num = one_message.begin[-1] + 1
    for i in range(0, one_message.n - len(one_message.begin)):
        one_message.begin.append(i_num)
        i_num += 1
    j_num = one_message.end[-1] + 1
    for i in range(0, one_message.n - len(one_message.end)):
        one_message.end.append(j_num)
        j_num += 1
    print("%s||%d||Go Cul, begin is %s, end is %s" %
        (time.strftime('%Y-%m-%d.%H:%M:%S',time.localtime(time.time())),comm_rank,one_message.begin,one_message.end))
    cycle_num = 0
    while not one_message.begin == one_message.end:
        cycle_num += 1
        #do somethings
        one_message.begin = update_begin_comb(one_message.begin, one_message.final)
    print("%s||%d||Cycle Number is %d" %
        (time.strftime('%Y-%m-%d.%H:%M:%S',time.localtime(time.time())),comm_rank,cycle_num))

def main():
    parser = _argparse()
    print("%s||%d||Initialization!" % (time.strftime('%Y-%m-%d.%H:%M:%S',time.localtime(time.time())),comm_rank))
    if comm_rank == 0:
        all_message = get_message(parser.n, parser.m, comm_size)
    else:
        all_message = None
    all_message = comm.scatter(all_message, root=0)
    print("%s||%d||begin: %s, end: %s, m: %d, n: %d" %
            (time.strftime('%Y-%m-%d.%H:%M:%S',time.localtime(time.time())), comm_rank, all_message.begin,
                all_message.end, all_message.m, all_message.n))
    #calculate function
    compound()
    print("%s||%d||Compelete END" % (time.strftime('%Y-%m-%d.%H:%M:%S',time.localtime(time.time())),comm_rank))

if __name__ == '__main__':
    main()
