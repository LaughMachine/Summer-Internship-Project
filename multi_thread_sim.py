import six_server_queue as ssq
import sys
import csv
import os
import time

from multiprocessing import Process

def writeLog(fil, table):
    c1 = csv.writer(fil)
    for val in table:
        c1.writerow(val)

arr = []
# Total_Time = 80000
# # Scale by nurses
Nurses = int(sys.argv[2])
# lbda_out = [1.0/(Nurses-Nurses*.1), 1.0/(Nurses-Nurses*.1)]
# mu_out = [1.0/2.0, 1.0/2.0]
# std_out = [1, 1]
# theta_out = [10000000, 10000000]
tau_out = [sys.argv[1], sys.argv[1]]
# k_out = 2
# hcost_out = [2,1]
# q_cap_out = [float('inf'), float('inf')]
# # Parallel simulation variables
# tot_par = 2
# s_alloc_out = [[Nurses/2,Nurses/2], [Nurses,Nurses]]
# rebalance1 = [1, 0]
# cont_out = [0, 1]
# preemption_out = [0, 1]
# time_vary = True

Total_Time = 80000
# Scale by nurses
# Nurses = 20
lbda_out = [1.0/(.24*Nurses), 1.0/(.24*Nurses)]
mu_out = [1.0/.5, 1.0/.5]
std_out = [1, 1]
theta_out = [10000, 10000]
# tau_out = [1, 1]
k_out = 2
hcost_out = [2,1]
q_cap_out = [float('inf'), float('inf')]
# Parallel simulation variables
tot_par = 3
s_alloc_out = [[Nurses/2,Nurses/2], [20,20], [Nurses/2, Nurses/2]]
rebalance1 = [1, 0, 0]
cont_out = [0, 1, 0]
preemption_out = [0, 1, 0]
time_vary = False
# Trial variables
trials = 10

s = ssq.Simulation(Total_Time, Nurses, lbda_out, mu_out, std_out, theta_out, tau_out, k_out, hcost_out, q_cap_out,
                   s_alloc_out, tot_par, rebalance1, cont_out, preemption_out, time_vary)
s.generate_arrivals(time_vary)
start_time = time.clock()
s.simulate(False,False)
end_time = time.clock()-start_time

out = []
out.append(s.arrival_count)
out.append(s.holding_cost)
out.append(s.time_server_occupied)
out.append([[x/s.t[ind] for ind, x in enumerate(y)] for y in s.weighted_ward])
out.append([[x/s.t[ind] for ind, x in enumerate(y)] for y in s.weighted_queue])
out.append([end_time])
print out
fil = open(os.getcwd() + "/Results" + str(sys.argv[3]) + ".csv", "wb")
writeLog(fil, out)