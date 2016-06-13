import numpy as np 
import heapq
import os
import time
import random 
import csv
import scipy as sp
import scipy.stats
import sys

class Patient:

    def __init__(self, t, t_arr, pt, location, aban, serv):
        self.t = t
        self.location = location
        self.ward = 0
        self.t_arr = t_arr
        self.pt = pt
        self.aban = aban
        self.serv = serv

    def __cmp__(self, other):
        return cmp(self.t, other.t)

    def get_time(self):
        return self.t

    def set_time(self, t):
        self.t = t

    def get_location(self):
        return self.location

    def set_location(self, location):
        self.location = location

    def get_pt(self):
        return self.pt

    def set_pt(self, pt):
        self.pt = pt

    def get_aban(self):
        return self.aban

    def get_serv(self):
        return self.serv 

    def set_serv(self, serv1):
        self.serv = serv1


class Simulation:

    def __init__(self, T, N, lbda, mu, std, theta, tau, classes, hcost, q_cap, s_alloc, par_sim, rb, cont, preemption, vary): 	
        self.pm = par_sim           # Total parallel simulations
        # ----------------- Environment Variables (same for all simulations) -----------------
        self.l_arr = lbda           # Arrival Rates
        self.l_aban = theta         # Abandonment Rates
        self.w_mu = mu              # Service Rate
        self.w_std = std            # STD of Service Rate
        self.k = classes            # Patient classes
        self.h_cost = hcost         # Holding cost for each patient type
        self.q_capac = q_cap        # Capacity of each queue
        self.N = N                  # Total Nurses
        self.Time = T               # Total Simulation Run-Time
        self.vary = vary            # Time varying arrival option
        self.rebal = rb             # Rebalance option
        # ----------------- Environment Variables (varies for each simulation) -----------------
        self.r_time = tau           # Shift Length
        self.preempt = preemption   # Preemption option
        self.cont = cont            # Continuous option
        # ----------------- Simulation Variables (arrays) -----------------
        self.arrival_list = []                              # Array containing the list of arrivals used for all simulations
        self.next_arrival = [0 for x in range(par_sim)]     # Index for of next arrival
        self.event_list = [[] for i in range(par_sim)]      # Arrays containing list of non-arrival events for each simulation
        self.t = [0 for i in range(par_sim)]                # Variable that tracks the current time of each simulation
        self.queue = [[[] for j in range(self.k)] for i in range(par_sim)]           # Empty arrays to hold patient objects in queue
        self.n_free = [N for i in range(par_sim)]           # Variable that tracks the current number of free nurses in each simulation
        self.r_time_arr = [i for i in self.r_time]               # Variable that tracks the next shift time change
        self.ward_alloc = [[j for j in s_alloc[i]] for i in range(par_sim)]         # Initiate arrays for the allocation of nurses to each ward
        self.ward_capac = [[j for j in s_alloc[i]] for i in range(par_sim)]         # Current Capacity of Each ward
        self.ward_nurse_def = [[0 for j in range(self.k)] for i in range(par_sim)]  # Keeps deficit of nurses during rebalance
        self.ward_assignment = [[] for i in range(par_sim)]                         # Keeps track of wards that still need nurses assigned during rebalance
        # ----------------- Simulation Variables (counters) -----------------
        self.balk_count = [0 for i in range(par_sim)]               # Initiate Balk count for each sim
        self.arrival_count = [0 for i in range(par_sim)]            # Initiate Arrival count for each sim
        self.abandonment_count = [0 for i in range(par_sim)]        # Initiate Abandonment count for each sim
        self.treated = [0 for i in range(par_sim)]                  # Initiate Treated patient count for each sim
        self.holding_cost = [0 for i in range(par_sim)]	            # Initiate holding cost
        self.time_server_occupied = [0 for i in range(par_sim)]     # Initiate time occupied
        # ----------------- Simulation Variables (patient class counters) -----------------
        self.weighted_ward = [[0 for j in range(self.k)] for i in range(par_sim)]   # Initate array for weighted ward count
        self.weighted_queue = [[0 for j in range(self.k)] for i in range(par_sim)]  # Initate array for weighted queue count
        self.queue_length = [[0 for j in range(self.k)] for i in range(par_sim)]    # Initiate array for queue length
        # ----------------- Other Variables -----------------
        self.statistics = [[] for i in range(par_sim)]

    # Method for randomly generating the list of arrivals
    def generate_arrivals(self, vary):
        arrivals, temp, t = [], [], 0
        # Generate initial set of arrivals 
        for pt, i in enumerate(self.l_arr):
            arr = vs_exp(1/float(i), 1, 1, 0, vary)
            temp.append(Patient(arr, 0, pt, 'arrival', np.random.exponential(self.l_aban[pt]), 
            np.random.exponential(self.w_mu[pt])))
        # Track a small set of arrivals and generate new ones as we push them to the overall arrival list
        heapq.heapify(temp)
        while (t < self.Time):
            arrivals.append(temp[0])
            arr = vs_exp(1/float(self.l_arr[pt]), 1, 1, 0, vary)
            t, pt = temp[0].get_time(), temp[0].get_pt()
            heapq.heappushpop(temp, Patient(arr + t, t, pt, 'arrival', np.random.exponential(self.l_aban[pt]), 
                np.random.exponential(self.w_mu[pt])))
        self.arrival_list = arrivals

    # Method for setting the list of arrivals to a predetermined list
    def set_arrivals(self, arrivals):
        self.arrival_list = arrivals

    # Method for randomly generating preexisting patients in the wards
    def set_preexisting_rnd(self, patient_types, serv_times, aban_times):
        existing = []
        for ind, e_p in enumerate(patient_types):
            existing.append(Patient(0, 0, e_p, 'arrival', aban_times[ind], serv_times[ind]))
        for i in range(self.pm):
            self.event_list[i] = [x for x in existing]

    # Method for running the simulation
    def simulate(self, use_preset_arr, save_data):
        for curr_sim in range(self.pm):
            while(self.t[curr_sim] < self.Time):
                hc, t_prev, servers_occupied, event_type,  = 0, 0, self.N - self.n_free[curr_sim], None
                ww, wq, ward_cnt, queue_cnt = [], [], [], []
                for wt in range(self.k):
                    hc += self.queue_length[curr_sim][wt]*self.h_cost[wt]
                    ward_cnt.append(self.ward_alloc[curr_sim][wt] - self.ward_capac[curr_sim][wt] + \
                                  self.queue_length[curr_sim][wt])
                    queue_cnt.append(self.queue_length[curr_sim][wt])
                    ww.append(ward_cnt[wt])
                    wq.append(queue_cnt[wt])
                # Get next event that will occur
                if self.event_list[curr_sim]:
                    curr_event_E = self.event_list[curr_sim][0]
                    curr_event_A = self.arrival_list[self.next_arrival[curr_sim]]
                    if curr_event_A.get_time() < curr_event_E.get_time():
                        curr_event = curr_event_A
                    else:
                        curr_event = curr_event_E
                else:
                    curr_event = self.arrival_list[self.next_arrival[curr_sim]]
                # Main body of the simulation
                t_prev = self.t[curr_sim]
                if self.cont[curr_sim] != 1 and self.rebal[curr_sim] == 1 and self.r_time_arr[curr_sim] < curr_event.get_time():
                    self.t[curr_sim] = self.r_time_arr[curr_sim]
                    self.r_time_arr[curr_sim] += self.r_time[curr_sim]
                    self._rebalance(curr_sim)
                    event_type = 'rebalance'
                else:
                    self.t[curr_sim] = curr_event.get_time()
                    if curr_event.get_location() == 'arrival':
                        self._arrival_event(curr_sim)
                        event_type = 'arrival'
                    elif curr_event.get_location() == 'ward':
                        self._departure_event(curr_sim)
                        event_type = 'departure'
                    elif curr_event.get_location() == 'abandonment':
                        self._aban_event(curr_sim)
                        event_type = 'abandonment'
                # Record length of queue and patients in wards and holding cost and servers occupied via a weighted sum
                self.holding_cost[curr_sim] += (self.t[curr_sim] - t_prev)*hc
                self.time_server_occupied[curr_sim] += (self.t[curr_sim] - t_prev)*servers_occupied
                for i in range(self.k):
                    self.weighted_queue[curr_sim][i] += (self.t[curr_sim] - t_prev)*queue_cnt[i]
                    self.weighted_ward[curr_sim][i] += (self.t[curr_sim] - t_prev)*ward_cnt[i]
                if save_data:
                    row = []
                    row.append(self.t[curr_sim])
                    row += ww + wq
                    row.append(event_type)
                    row.append([x for x in self.ward_alloc[curr_sim]])
                    row.append([x for x in self.ward_capac[curr_sim]])
                    row.append(self.n_free[curr_sim])
                    row.append(len(self.event_list[curr_sim]))
                    self.statistics[curr_sim].append(row)

    def _arrival_event(self, sim):
        self.arrival_count[sim] += 1
        pt = self.arrival_list[self.next_arrival[sim]].get_pt()
        serv_time = self.arrival_list[self.next_arrival[sim]].get_serv()
        aban_time = self.arrival_list[self.next_arrival[sim]].get_aban()
        self.next_arrival[sim] += 1
        admitted = False
        # Check if patient can be admitted in case where nurses are free
        if self.n_free[sim] > 0:
            admitted = True if self.cont[sim] == 1 else True if self.ward_capac[sim][pt] > 0 else False
        # Check if there is preemption in the case patient hasn't been admited
        if not admitted:
            if self.preempt[sim] == 1 and self.cont[sim] == 1:
                # Find patient in ward with smallest cost
                min_cost = self.h_cost[pt]
                min_pt = pt
                for h in range(self.k):
                    if self.h_cost[h] < min_cost and (self.ward_alloc[sim][h]-self.ward_capac[sim][h]) > 0:
                        min_cost, min_pt = self.h_cost[h], h
                # If patient in ward with smaller cost is found, remove push back into queue
                if self.h_cost[pt] > min_cost:
                    for pt_ind, pt_var in reversed(list(enumerate(self.event_list[sim]))):
                        if pt_var.get_pt() == min_pt:
                            pt0 = pt_var
                            del self.event_list[sim][pt_ind]
                            heapq.heapify(self.event_list[sim])
                            break
                    pt0.set_serv(pt0.get_time()-self.t[sim])
                    pt0.set_location('abandonment')
                    self.queue[sim][min_pt].insert(0,pt0)
                    self.queue_length[sim][min_pt] += 1
                    self.ward_capac[sim][min_pt] += 1
                    self.n_free[sim] += 1
                    admitted = True
        if admitted:
            new_ward_arr = Patient(serv_time + self.t[sim], self.t[sim], pt, 'ward', 
                aban_time, serv_time)
            heapq.heappush(self.event_list[sim], new_ward_arr)
            self.ward_capac[sim][pt] -= 1
            self.n_free[sim] -= 1 
        else:
            if self.queue_length[sim][pt] >= self.q_capac[pt]:
                self.balk_count += 1
            else:
                new_queue_arr = Patient(aban_time + self.t[sim], self.t[sim], pt, 'abandonment', 
                    aban_time, serv_time)
                self.queue[sim][pt].append(new_queue_arr)
                self.queue_length[sim][pt] += 1

    def _departure_event(self, sim):
        # Increase capacity in ward where patient is leaving
        pt = self.event_list[sim][0].get_pt()
        self.ward_capac[sim][pt] += 1
        heapq.heappop(self.event_list[sim])
        self.treated[sim] += 1
        new_patient = False
        # Find which patient in queue will replace the outgoing patient (get correct pt)
        if self.cont[sim] == 1:
            max_cost = 0
            # Find the max cost patient in waiting
            for h in range(self.k):
                if self.h_cost[h] >= max_cost and self.queue_length[sim][h] > 0:
                    max_cost = self.h_cost[h]
                    pt = h
            # check initial patient class has a queue
            if self.queue_length[sim][pt] > 0:
                new_patient = True
            else:
                self.n_free[sim] += 1
        else:
            # Check if patients are in queue, checks if ward capacity allows for another patient
            if self.queue[sim][pt] and self.ward_capac[sim][pt] > 0 and self.n_free[sim] >= 0:
                new_patient = True
            # Case when nurse finishes treating a patient and moves to newly assigned ward
            elif self.ward_capac[sim][pt] <= 0 and self.n_free[sim] >= 0 and self.ward_assignment[sim]:
                pt = random.choice(self.ward_assignment[sim])
                self.ward_nurse_def[sim][pt] += 1
                if self.ward_nurse_def[sim][pt] == 0:
                    self.ward_assignment[sim].remove(pt)
                if self.queue_length[sim][pt] > 0:
                    new_patient = True
                else:
                    self.n_free[sim] += 1
            else:
                self.n_free[sim] += 1
        # If there is a new patient push to event_list and adjust ward status
        if new_patient:
            next_patient = self.queue[sim][pt].pop(0)
            self.queue_length[sim][pt] -= 1 
            self.ward_capac[sim][pt] -= 1
            next_patient.set_time(self.t[sim] + next_patient.get_serv())
            next_patient.set_location('ward')
            heapq.heappush(self.event_list[sim], next_patient)

    def _rebalance(self, sim):
        old_alloc = [x for x in self.ward_alloc[sim]]
        self._set_new_alloc(sim, old_alloc)
        if self.preempt[sim] == 1:
            self._reset_wards(sim)
        self._fill(sim)
        for i in range(self.k):
            if self.n_free[sim] == 0 and self.ward_capac[sim][i] > 0:
                self.ward_nurse_def[sim][i] = -self.ward_capac[sim][i]
                self.ward_assignment[sim].append(i)

    def _fill(self, sim):
        # For each type of class of patients we will check to see if we can refill the wards
        for i in range(self.k):
            # If there are nurses free and the ward has capacity and has a queue, we assign patients
            while(self.n_free[sim] > 0 and self.ward_capac[sim][i] > 0 and self.queue_length[sim][i] > 0):
                # Get next patient in the queue
                next_patient = self.queue[sim][i].pop(0)
                self.queue_length[sim][i] -= 1
                # Change patient type to ward patient from abandoner
                next_patient.set_time(self.t[sim] + next_patient.get_serv())
                next_patient.set_location('ward')
                # Add patient to Events list
                heapq.heappush(self.event_list[sim], next_patient)
                # Update counters 
                self.n_free[sim] -= 1
                self.ward_capac[sim][i] -= 1

    def _reset_wards(self, sim):
        requeue = []
        for i in range(self.k):
            requeue.append([])
        # Changes the ward patients back to queue patients after readjusting service time
        # to remaining service time. 
        for pat in self.event_list[sim]:
            remaining_service = pat.get_time() - self.t[sim]
            pat.set_serv(remaining_service)
            pat.set_location('abandonment')
            requeue[pat.get_pt()].append(pat)
        # Reset the wards so they are empty and nurses are all freed
        self.event_list[sim] = []
        self.ward_capac[sim] = [x for x in self.ward_alloc[sim]]
        self.n_free[sim] = self.N
        # Procedure pushes all patients in wards back to their respective queues
        for i in range(self.k):
            incoming_length = len(requeue[i])
            self.queue[sim][i] = requeue[i] + self.queue[sim][i]
            self.queue_length[sim][i] += incoming_length

    def _aban_event(self, sim):
        pt = self.event_list[sim][0].get_pt()
        abandoner = heapq.heappop(self.event_list[sim])
        if abandoner in self.queue[sim][pt]:
            self.queue[sim][pt].remove(abandoner)
            self.abandonment_count[sim] += 1
            self.queue_length[sim][pt] -= 1

    def _set_new_alloc(self, sim, old_alloc):
        s, k, N = 0, self.k, self.N
        total = max(sum(self.queue_length[sim]), k)
        for i in range(k-1):
            old_alloc[i] = self.ward_alloc[sim][i]
            self.ward_alloc[sim][i] = self.queue_length[sim][i]*(N-k)/total + 1
            self.ward_capac[sim][i] = self.ward_capac[sim][i] + self.ward_alloc[sim][i] - old_alloc[i]
            s += self.queue_length[sim][i]*(N-k)/total
        # Setting new allocation and calculating new capacities
        old_alloc[k-1] = self.ward_alloc[sim][k-1]
        self.ward_alloc[sim][k-1] = (N - k) - s + 1
        self.ward_capac[sim][k-1] = self.ward_capac[sim][k-1] + self.ward_alloc[sim][k-1] - old_alloc[k-1]

# ------------------- Numerical Methods -------------------
def mean_confidence_interval(data, confidence=0.95): 
    a = 1.0*np.array(data)
    n = len(a)
    m, se = np.mean(a), scipy.stats.sem(a)
    h = se * sp.stats.t._ppf((1+confidence)/2., n-1)
    return m, m-h, m+h

def vs_exp(l, period, a, t, vary):
    if vary:
        original_t = t
        lv = 2*a*float(l)
        U = 1.1
        while(U > lv/(2*a*float(l))):
            U = np.random.uniform()
            t = t + np.random.exponential(1/float(l*a*2))
            lv = l*a + a*l*np.sin(2*np.pi*t/period)
        return (t - original_t)
    else:
        return np.random.exponential(1/l)

def writeLog(fil, table):
    c1 = csv.writer(fil)
    for val in table:
        c1.writerow(val)

if __name__ == "__main__":
    # ================ Input Variables ================
    Total_Time = 400
    lbda_out = [1.0/15.0, 1.0/15.0]
    mu_out = [1.0/8.0, 1.0/8.0]
    std_out = [1, 1]
    theta_out = [10000, 10000]
    tau_out = [.5, .5]
    k_out = 2
    Nurses = 4
    hcost_out = [1,2]
    q_cap_out = [float('inf'), float('inf')]

    # Parallel simulation variables
    s_alloc_out = [[2,2], [4,4]]
    rebalance1 = [1, 0]
    cont_out = [0, 1]
    preemption_out = [0, 1]
    time_vary = [False, False]

    s = Simulation(Total_Time, Nurses, lbda_out, mu_out, std_out, theta_out, tau_out, k_out, hcost_out, q_cap_out,
                   s_alloc_out, 2, rebalance1, cont_out, preemption_out, time_vary)
    s.generate_arrivals(False)
    s.simulate(False, True)
    print s.arrival_count
    print len(s.arrival_list)
    print s.holding_cost
    print s.time_server_occupied
    print [[y/s.t[ind] for y in x] for ind, x in enumerate(s.weighted_queue)]
    print [[y/s.t[ind] for y in x] for ind, x in enumerate(s.weighted_ward)]

    fil0 = open(os.getcwd() + "/Sim_Rebalance_6.csv", "wb")
    fil1 = open(os.getcwd() + "/Sim_No_Rebalance_6.csv", "wb")

    writeLog(fil0, s.statistics[0])
    writeLog(fil1, s.statistics[1])

