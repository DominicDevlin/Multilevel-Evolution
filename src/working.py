import os
import sys
import csv
import math

# what type of simulation, 0 = V, 1 = m * V, 2 = mV * inits
sim_type = 2

def run_simulation2(threshold, mutation, init):
    kp = init[0]
    dp = init[1]
    kq = init[2]
    dq = init[3]
    name = threshold + "-" + mutation + "-" + kp + dp + kq + dq
    os.system("g++ -O3 fission.cpp -o " + name + " -fopenmp")
    os.system("./" + name + " " + threshold + " " + mutation + " " + kp + " " + dp + " " + kq + " " + dq)
    
def run_simulation1(threshold, mutation):
    name = threshold + "-" + mutation
    os.system("g++ -O3 fission.cpp -o " + name + " -fopenmp")
    os.system("./" + name + " " + threshold + " " + mutation)
    
def run_simulation0(threshold):
    os.system("g++ -O3 fission.cpp -o " + threshold + " -fopenmp")
    os.system("./" + threshold + " " + threshold)
    

thresholds = []
mutations = []
n_thresholds = 0
n_inits = 0

inits = []


     
with open('params.csv', mode = 'r') as f:
    csvFile = csv.reader(f) 
    for line in csvFile:
        if (sim_type > 0):
            thresholds.append(line[0])
            mutations.append(line[1])
            n_thresholds += 1
        else:
            thresholds.append(line[0])
            n_thresholds += 1
            

if (sim_type == 2):
    with open('init.csv', mode = 'r') as f:
        csvFile = csv.reader(f) 
        for line in csvFile:
            inits.append(line)
            n_inits += 1


if (sim_type == 2):
    index = int(sys.argv[1])
    ind1 = math.floor(index / n_inits)
    ind2 = index % n_inits
    thresh = thresholds[ind1]
    mut = mutations[ind1]
    init = inits[ind2]
    run_simulation2(thresh, mut, init)
elif (sim_type == 1):
    index = int(sys.argv[1])
    thresh = thresholds[index]
    mut = mutations[index]
    run_simulation1(thresh, mut)
else:
    index = int(sys.argv[1])
    thresh = thresholds[index]
    run_simulation0(thresh)








