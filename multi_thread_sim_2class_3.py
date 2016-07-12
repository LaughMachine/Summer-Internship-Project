#!/apps/anaconda-2.3.0/bin/python

import six_server_queue as ssq
import sys
import csv
import os
import time

def writeLog(fil, table):
    c1 = csv.writer(fil)
    for val in table:
        c1.writerow(val)

raw_data = []
with open(sys.argv[1], 'rb') as fil:
    r = csv.reader(fil)
    for row in r:
        raw_data.append(row)

arr = []
ratio = [1, 3]
tot_par = int(raw_data[3][1])
Nurses = int(sys.argv[3])
tau_out = [int(sys.argv[2]) for x in range(tot_par)]
k_out = int(raw_data[8][1])
Total_Time = int(raw_data[0][1])
lbda_out = [1.0/(float(sys.argv[4])*float(raw_data[9][x])*Nurses*ratio[x-1]/float(sum(ratio))) for x in range(1,k_out+1)]
mu_out = [1.0/float(raw_data[9][x]) for x in range(1,k_out+1)]
std_out = [float(raw_data[10][x]) for x in range(1,k_out+1)]
theta_out = [float(raw_data[11][x]) for x in range(1,k_out+1)]
hcost_out = [float(raw_data[12][x]) for x in range(1,k_out+1)]
q_cap_out = [float(raw_data[13][x]) for x in range(1,k_out+1)]
rebalance1 = [int(raw_data[4][x]) for x in range(1,tot_par+1)]
cont_out = [int(raw_data[5][x]) for x in range(1,tot_par+1)]
preemption_out = [int(raw_data[6][x]) for x in range(1,tot_par+1)]
time_vary = False if int(raw_data[7][1]) == 0 else True
s_alloc_out = []
dedicated_alloc = [int(Nurses*.27), Nurses-int(Nurses*.27)]
for i in range(0, tot_par):
    if cont_out[i] == 1:
        s_alloc_out.append([Nurses, Nurses])
    else:
        s_alloc_out.append(dedicated_alloc)


print 'Nurses: ' + str(Nurses)
print 'tau_out: ' + str(tau_out)
print 'classes: ' + str(k_out)
print 'total time: ' + str(Total_Time)
print 'arrival rate: ' + str(lbda_out)
print 'service rate: ' + str(mu_out)
print 'service std: ' + str(std_out)
print 'theta: ' + str(theta_out)
print 'holding cost: ' + str(hcost_out)
print 'queue length: ' + str(q_cap_out)
print 'parallel simulations' + str(tot_par)
print 'starting allocation: ' + str(s_alloc_out)
print 'rebalance: ' + str(rebalance1)
print 'continuous: ' + str(cont_out)
print 'preemption: ' + str(preemption_out)
print 'time vary: ' + str(time_vary)


s = ssq.Simulation(Total_Time, Nurses, lbda_out, mu_out, std_out, theta_out, tau_out, k_out, hcost_out, q_cap_out,
                   s_alloc_out, tot_par, rebalance1, cont_out, preemption_out, time_vary)
s.generate_arrivals(time_vary)
start_time = time.clock()
s.simulate(False,False)
end_time = time.clock()-start_time

out = []
out.append([x/float(s.t[ind]) for ind, x in enumerate(s.arrival_count)])
out.append([x/float(s.t[ind]) for ind, x in enumerate(s.holding_cost)])
out.append([x/float(Nurses*s.t[ind]) for ind, x in enumerate(s.time_server_occupied)])
out.append([[x/float(s.t[ind]) for ind, x in enumerate(y)] for y in s.weighted_ward])
out.append([[x/float(s.t[ind]) for ind, x in enumerate(y)] for y in s.weighted_queue])
out.append([end_time])


fil = open(os.getcwd() + "/Results" + str(sys.argv[5]) + ".csv", "wb")
writeLog(fil, out)