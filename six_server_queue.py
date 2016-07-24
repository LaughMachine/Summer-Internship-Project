#!/apps/anaconda-2.3.0/bin/python

import numpy as np
import heapq
import os
import time
import random 
import csv
import scipy as sp
import scipy.stats
from scipy.integrate import odeint
from scipy.optimize import brentq
from scipy.optimize import fmin_tnc
from scipy.optimize import minimize
from math import factorial

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
        self.safety = 1             # Safety variable
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
        self.dedicated_alloc = [find_dedicated_2_class(self.N,self.l_arr,self.w_mu,
                                            self.h_cost) for i in range(par_sim)]   # Keep track of dedicated capacity allocation
        self.ward_alloc = [[j for j in self.dedicated_alloc[i]] for i in
                           range(par_sim)]  # Initiate arrays for the allocation of nurses to each ward
        self.ward_capac = [[j for j in self.dedicated_alloc[i]] for i in range(par_sim)]  # Current Capacity of Each ward
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
        self.alloc_at_rebal = [[] for i in range(par_sim)]
        self.head_count_at_rebal = [[] for i in range(par_sim)]
        self.dedicated_cost = [cost_func(self.dedicated_alloc[i],self.l_arr,self.w_mu,self.N,self.h_cost) for i in range(par_sim)]

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
            t, pt = temp[0].get_time(), temp[0].get_pt()
            arr = vs_exp(1/float(self.l_arr[pt]), 1, 1, t, vary)
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
                    # Number of patients in the ward and queue for a specific class
                    ward_cnt.append(self.ward_alloc[curr_sim][wt] - self.ward_capac[curr_sim][wt] + \
                                  self.queue_length[curr_sim][wt])
                    # Number of patients in the queue for a specific class
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
                if self.cont[curr_sim] != 1 and self.rebal[curr_sim] >= 1 and self.r_time_arr[curr_sim] < curr_event.get_time():
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
            self.event_list[curr_sim] = []
                # print str(curr_sim) + ' ' + str(self.t[curr_sim])

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
                min_cost = self.h_cost[pt]/float(self.w_mu[pt])
                min_pt = pt
                for h in range(self.k):
                    if self.h_cost[h]/float(self.w_mu[h]) < min_cost and (self.ward_alloc[sim][h]-self.ward_capac[sim][h]) > 0:
                        min_cost, min_pt = self.h_cost[h]/float(self.w_mu[h]), h
                # If patient in ward with smaller cost is found, remove push back into queue
                if self.h_cost[pt]/float(self.w_mu[pt]) > min_cost:
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
                if self.h_cost[h]/float(self.w_mu[h]) >= max_cost and self.queue_length[sim][h] > 0:
                    max_cost = self.h_cost[h]/float(self.w_mu[h])
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
        if self.rebal[sim] == 1:
            new_alloc = self._get_new_alloc_ode(sim, old_alloc)
        elif self.rebal[sim] == 2:
            new_alloc = self._get_new_alloc_gen(sim, old_alloc)
        elif self.rebal[sim] == 3:
            norm_lbda = [1 / x for x in self.l_arr]
            mu = [1 / x for x in self.w_mu]
            sum = 0
            for l in range(self.k):
                sum += norm_lbda[l]/(float(self.N)*mu[l])
            if sum <.94:
                new_alloc = self._get_new_alloc_multi_heur_92(sim)
            else:
                new_alloc = self._get_new_alloc_multi_heur_96(sim)
        elif self.rebal[sim] == 4:
            new_alloc = self._get_new_alloc_ode_d1(sim, self.safety,0)
        elif self.rebal[sim] == 5:
            new_alloc = self._get_new_alloc_ode_d2(sim, self.safety,0)
        elif self.rebal[sim] == 6:
            new_alloc = self._get_new_alloc_ode_d1(sim, self.safety, 1)
        elif self.rebal[sim] == 7:
            new_alloc = self._get_new_alloc_ode_d2(sim, self.safety, 1)
        elif self.rebal[sim] == 8:
            new_alloc = self._get_new_alloc_three_period_d1_tnc(sim, self.safety)
        elif self.rebal[sim] == 9:
            new_alloc = self._get_new_alloc_three_period_d2_tnc(sim, self.safety)
        elif self.rebal[sim] == 10:
            new_alloc = self._get_new_alloc_three_period_d3_tnc(sim, self.safety)
        elif self.rebal[sim] == 11:
            new_alloc = self._get_new_alloc_three_period_d1_slsqp(sim, self.safety, 0)
        elif self.rebal[sim] == 12:
            new_alloc = self._get_new_alloc_three_period_d2_slsqp(sim, self.safety, 0)
        elif self.rebal[sim] == 13:
            new_alloc = self._get_new_alloc_three_period_d1_slsqp(sim, self.safety, 1)
        elif self.rebal[sim] == 14:
            new_alloc = self._get_new_alloc_three_period_d2_slsqp(sim, self.safety, 1)
        else:
            print 'error no rebalance policy'
            new_alloc = old_alloc
        self.alloc_at_rebal[sim].append(new_alloc)
        head_count = [self.ward_alloc[sim][i] - self.ward_capac[sim][i] + self.queue_length[sim][i] for i in range(self.k)]
        self.head_count_at_rebal[sim].append(head_count)
        # Setting new allocation and calculating new capacities
        for i in range(self.k):
            self.ward_alloc[sim][i] = new_alloc[i]
            self.ward_capac[sim][i] = self.ward_capac[sim][i] + self.ward_alloc[sim][i] - old_alloc[i]

    def _get_new_alloc_gen(self, sim, old_alloc):
        s, k, N = 0, self.k, self.N
        new_alloc = []
        total = max(sum(self.queue_length[sim]), k)
        for i in range(k-1):
            new_alloc.append(self.queue_length[sim][i]*(N-k)/total + 1)
            s += self.queue_length[sim][i] * (N - k) / total
        new_alloc.append((N - k) - s + 1)
        return new_alloc

    def _get_new_alloc_ode_d2(self, sim, safe_coef, rnd):
        y0 = []
        t0 = [0 for x in range(self.k)]
        new_alloc = []
        norm_lbda = [1/(x*float(self.N)) for x in self.l_arr]
        mu = [1/x for x in self.w_mu]
        safety = [safe_coef*self.N**.5, 0]
        for i in range(self.k):
            y0.append(max((self.ward_alloc[sim][i] - self.ward_capac[sim][i] + self.queue_length[sim][i] - safety[i])/float(self.N),0))
        if sum(y0) > 1:
            if y0[0] >= (mu[0]-norm_lbda[0])*self.r_time[0]+1:
                return [self.N, 0]
            else:
                solution = ode_sys_complete(self.r_time[sim], y0, t0, norm_lbda, mu, self.k)
                if rnd == 0:
                    for i in solution:
                        new_alloc.append(int(np.around(self.N*i/self.r_time[sim])))
                    new_alloc = rounding_pref_2_class(new_alloc, self.N)
                    return new_alloc
                elif rnd == 1:
                    u1 = np.around(solution[0] * self.N / self.r_time[sim])
                    return [u1, self.N - u1]
                else:
                    print 'error'
        else:
            return self.dedicated_alloc[sim]

    def _get_new_alloc_ode_d1(self, sim, safe_coef, rnd):
        y0 = []
        t0 = [0 for x in range(self.k)]
        new_alloc = []
        norm_lbda = [1 / (x * float(self.N)) for x in self.l_arr]
        mu = [1 / x for x in self.w_mu]
        rho = [norm_lbda[i] / mu[i] for i in range(self.k)]
        safety = [safe_coef*self.N**.5, 0]
        for i in range(self.k):
            y0.append(max(
                (self.ward_alloc[sim][i] - self.ward_capac[sim][i] + self.queue_length[sim][i] - safety[i]) / float(
                    self.N), 0))
        if y0[0] <= rho[0] and y0[1] <= rho[1]:
            return self.dedicated_alloc[sim]
        else:
            if y0[0] >= (mu[0] - norm_lbda[0]) * self.r_time[0] + 1:
                return [self.N, 0]
            else:
                solution = ode_sys_complete(self.r_time[sim], y0, t0, norm_lbda, mu, self.k)
                if rnd == 0:
                    for i in solution:
                        new_alloc.append(int(np.around(self.N * i / self.r_time[sim])))
                    new_alloc = rounding_pref_2_class(new_alloc, self.N)
                    return new_alloc
                elif rnd == 1:
                    u1 = np.around(solution[0] * self.N / self.r_time[sim])
                    return [u1, self.N - u1]
                else:
                    print 'error'

    def _get_new_alloc_ode(self, sim, old_alloc):
        y0 = []
        t0 = [0 for x in range(self.k)]
        new_alloc = []
        norm_lbda = [1 / (x * float(self.N)) for x in self.l_arr]
        mu = [1 / x for x in self.w_mu]
        rho = [norm_lbda[i] / mu[i] for i in range(self.k)]
        safety = [self.N ** .5, 0]
        # safety = [0, 0]
        for i in range(self.k):
            y0.append(max(
                (self.ward_alloc[sim][i] - self.ward_capac[sim][i] + self.queue_length[sim][i] - safety[i]) / float(
                    self.N), 0))
        if y0[0] <= rho[0] and y0[1] <= rho[1]:
            return self.dedicated_alloc[sim]
        else:
        # if sum(y0) > 1:
            if y0[0] >= (mu[0] - norm_lbda[0]) * self.r_time[0] + 1:
                return [self.N, 0]
            else:
                solution = ode_sys_complete(self.r_time[sim], y0, t0, norm_lbda, mu, self.k)
                for i in solution:
                    new_alloc.append(int(np.around(self.N * i / self.r_time[sim])))
                new_alloc = rounding_pref_2_class(new_alloc, self.N)
                # u1 = np.around(solution[0]*self.N/self.r_time[sim])

                # while sum(new_alloc) != self.N:
                #     if sum(new_alloc) > self.N:
                #         if new_alloc[1] > 0:
                #             new_alloc[1] -= 1
                #         else:
                #             new_alloc[0] -= 1
                #     else:
                #         new_alloc[0] += 1
                return new_alloc
                # return [u1, self.N - u1]
        # else:
        #     return self.dedicated_alloc[sim]

    def _get_new_alloc_multi_heur_0(self, sim):
        safety = 2 * self.N**.5
        x0 = [(self.ward_alloc[sim][i] - self.ward_capac[sim][i] + self.queue_length[sim][i])/float(self.N)
              for i in range(self.k)]
        norm_lbda = [1 / (x * float(self.N)) for x in self.l_arr]
        mu = [1 / x for x in self.w_mu]
        rho = [norm_lbda[i] / mu[i] for i in range(self.k)]
        # Both are below utilization
        if x0[0] <= rho[0] and x0[1] <= rho[1]:
            return self.dedicated_alloc[sim]
        # Case where atleast one is above utilization
        else:
            u_bar = [(x0[i] + norm_lbda[i] * self.r_time[sim]) / (1 + mu[i] * self.r_time[sim]) for i in range(self.k)]
            new_u_bar = [(x0[i] - safety + norm_lbda[i] * self.r_time[sim]) / (1 + mu[i] * self.r_time[sim]) for i
                         in range(self.k)]
            # Case where ward 1 is over capacity and ward 2 is not
            if x0[0] > rho[0] and x0[1] <= rho[1]:
                u1 = max(self.dedicated_alloc[sim][0], int(np.around(self.N * min(1, new_u_bar[0]))))
                return [u1, self.N - u1]
            # Case where ward 2 is over capacity and ward 1 is not
            elif x0[0] <= rho[0] and x0[1] > rho[1]:
                u1 = max(np.ceil(self.N * rho[0]), self.dedicated_alloc[sim][0])
                return [u1, self.N - u1]
            # Case where both ward 1 and 2 are over capacity
            else:
                # Case where we can empty both 1 and 2
                if sum(u_bar) <= 1:
                    root = find_foc_root(x0, u_bar, self.h_cost, self.r_time[sim], mu, norm_lbda)
                    u1 = int(np.around(root * self.N))
                    return [u1, self.N - u1]
                # Case where we cannot empty both, try to drain 1
                else:
                    u1 = max(self.dedicated_alloc[sim][0], int(np.around(self.N * min(1, new_u_bar[0]))))
                    return [u1, self.N - u1]

    def _get_new_alloc_multi_heur_1(self, sim):
        safety = (self.N ** .5)
        x0 = [(self.ward_alloc[sim][i] - self.ward_capac[sim][i] + self.queue_length[sim][i] - safety) / float(self.N)
              for i in range(self.k)]
        norm_lbda = [1 / (x * float(self.N)) for x in self.l_arr]
        mu = [1 / x for x in self.w_mu]
        rho = [norm_lbda[i] / mu[i] for i in range(self.k)]
        # Both are below utilization
        if x0[0] <= rho[0] and x0[1] <= rho[1]:
            return self.dedicated_alloc[sim]
        # Case where atleast one is above utilization
        else:
            u_bar = [(x0[i] + norm_lbda[i] * self.r_time[sim]) / (1 + mu[i] * self.r_time[sim]) for i in range(self.k)]
            new_u_bar = [(x0[i] - safety + norm_lbda[i] * self.r_time[sim]) / (1 + mu[i] * self.r_time[sim]) for i
                         in range(self.k)]
            # Case where ward 1 is over capacity and ward 2 is not
            if x0[0] > rho[0] and x0[1] <= rho[1]:
                u1 = max(self.dedicated_alloc[sim][0], int(np.around(self.N * min(1, u_bar[0]))))
                return [u1, self.N - u1]
            # Case where ward 2 is over capacity and ward 1 is not
            elif x0[0] <= rho[0] and x0[1] > rho[1]:
                u1 = max(np.ceil(self.N * rho[0]), self.dedicated_alloc[sim][0])
                return [u1, self.N - u1]
            # Case where both ward 1 and 2 are over capacity
            else:
                # Case where we can empty both 1 and 2
                if sum(u_bar) <= 1:
                    root = find_foc_root(x0, u_bar, self.h_cost, self.r_time[sim], mu, norm_lbda)
                    u1 = max(int(np.around(root * self.N)), self.dedicated_alloc[sim][0])
                    # return [u1, self.N - u1]
                    return self.dedicated_alloc[sim]
                # Case where we cannot empty both, try to drain 1
                else:
                    u1 = max(self.dedicated_alloc[sim][0], int(np.around(self.N * min(1, u_bar[0]))))
                    return [u1, self.N - u1]

    def _get_new_alloc_multi_heur_92(self, sim):
        mu = [1 / x for x in self.w_mu]
        safety = [float((1 + mu[0] * self.r_time[sim]) * (self.N ** 0.5) * (2 / 3)) / float(self.N), 0]
        x0 = [(self.ward_alloc[sim][i] - self.ward_capac[sim][i] + self.queue_length[sim][i]) / float(self.N)
              for i in range(self.k)]
        norm_lbda = [1 / (x * float(self.N)) for x in self.l_arr]
        rho = [norm_lbda[i] / mu[i] for i in range(self.k)]
        # Both are below utilization
        if x0[0] <= rho[0] and x0[1] <= rho[1]:
            return self.dedicated_alloc[sim]
        # Case where atleast one is above utilization
        else:
            u_bar = [(x0[i] + norm_lbda[i] * self.r_time[sim]) / (1 + mu[i] * self.r_time[sim]) for i in range(self.k)]
            new_u_bar = [(x0[i] - safety[i] + norm_lbda[i] * self.r_time[sim]) / (1 + mu[i] * self.r_time[sim]) for i
                         in range(self.k)]
            # Case where ward 1 is over capacity and ward 2 is not
            if x0[0] > rho[0] and x0[1] <= rho[1]:
                u1 = max(self.dedicated_alloc[sim][0], int(np.around(self.N * min(1, new_u_bar[0]))))
                return [u1, self.N - u1]
            # Case where ward 2 is over capacity and ward 1 is not
            elif x0[0] <= rho[0] and x0[1] > rho[1]:
                # u1 = max(np.ceil(self.N * rho[0]), self.dedicated_alloc[sim][0])
                root = find_foc_root(x0, u_bar, self.h_cost, self.r_time[sim], mu, norm_lbda)
                u1 = int(np.around(root * self.N))
                return [u1, self.N - u1]
            # Case where both ward 1 and 2 are over capacity
            else:
                # Case where we can empty both 1 and 2
                if sum(u_bar) <= 1:
                    root = find_foc_root(x0, u_bar, self.h_cost, self.r_time[sim], mu, norm_lbda)
                    # u1 = int(np.around(root * self.N))
                    u1 = max(int(np.around(root * self.N)), self.dedicated_alloc[sim][0])
                    return [u1, self.N - u1]
                # Case where we cannot empty both, try to drain 1
                else:
                    u1 = max(self.dedicated_alloc[sim][0], int(np.around(self.N * min(1, new_u_bar[0]))))
                    return [u1, self.N - u1]

    def _get_new_alloc_multi_heur_96(self, sim):
        mu = [1 / x for x in self.w_mu]
        safety = [(1 + mu[0] * self.r_time[sim]) * (self.N ** 0.5) / float(self.N), 0]
        x0 = [(self.ward_alloc[sim][i] - self.ward_capac[sim][i] + self.queue_length[sim][i]) / float(self.N)
              for i in range(self.k)]
        norm_lbda = [1 / (x * float(self.N)) for x in self.l_arr]
        rho = [norm_lbda[i] / mu[i] for i in range(self.k)]
        # Both are below utilization
        if x0[0] <= rho[0] and x0[1] <= rho[1]:
            return self.dedicated_alloc[sim]
        # Case where atleast one is above utilization
        else:
            u_bar = [(x0[i] + norm_lbda[i] * self.r_time[sim]) / (1 + mu[i] * self.r_time[sim]) for i in range(self.k)]
            new_u_bar = [(x0[i] - safety[i] + norm_lbda[i] * self.r_time[sim]) / (1 + mu[i] * self.r_time[sim]) for i
                         in range(self.k)]
            # Case where ward 1 is over capacity and ward 2 is not
            if x0[0] > rho[0] and x0[1] <= rho[1]:
                u1 = max(self.dedicated_alloc[sim][0], int(np.around(self.N * min(1, new_u_bar[0]))))
                return [u1, self.N - u1]
            # Case where ward 2 is over capacity and ward 1 is not
            elif x0[0] <= rho[0] and x0[1] > rho[1]:
                # u1 = max(np.ceil(self.N * rho[0]), self.dedicated_alloc[sim][0])
                root = find_foc_root(x0, u_bar, self.h_cost, self.r_time[sim], mu, norm_lbda)
                u1 = int(np.around(root * self.N))
                return [u1, self.N - u1]
            # Case where both ward 1 and 2 are over capacity
            else:
                # Case where we can empty both 1 and 2
                if sum(u_bar) <= 1:
                    root = find_foc_root(x0, u_bar, self.h_cost, self.r_time[sim], mu, norm_lbda)
                    # u1 = int(np.around(root * self.N))
                    u1 = max(int(np.around(root * self.N)), self.dedicated_alloc[sim][0])
                    return [u1, self.N - u1]
                # Case where we cannot empty both, try to drain 1
                else:
                    u1 = max(self.dedicated_alloc[sim][0], int(np.around(self.N * min(1, new_u_bar[0]))))
                    return [u1, self.N - u1]

    def _get_new_alloc_single_period(self, sim):
        x0 = [(self.ward_alloc[sim][i] - self.ward_capac[sim][i] + self.queue_length[sim][i]) / float(self.N)
              for i in range(self.k)]
        norm_lbda = [1 / (x * float(self.N)) for x in self.l_arr]
        mu = [1 / x for x in self.w_mu]
        rho = [norm_lbda[i] / mu[i] for i in range(self.k)]

    def _get_new_alloc_three_period_d1_tnc(self, sim, safe_coef):
        mu = [1 / x for x in self.w_mu]
        safety = [safe_coef*self.N**.5, 0]
        x0 = [max((self.ward_alloc[sim][i] - self.ward_capac[sim][i] + self.queue_length[sim][i] - safety[i]) / float(self.N),0)
              for i in range(self.k)]
        norm_lbda = [1 / (x * float(self.N)) for x in self.l_arr]
        rho = [norm_lbda[i] / mu[i] for i in range(self.k)]
        u_bar = [(x0[i] + norm_lbda[i] * self.r_time[sim]) / (1 + mu[i] * self.r_time[sim]) for i in range(self.k)]
        if x0[0] <= rho[0] and x0[1] <= rho[1]:
            return self.dedicated_alloc[sim]
        elif x0[0] >= (mu[0]-norm_lbda[0])*self.r_time[0]+1:
            return [self.N, 0]
        else:
            new_result = fmin_tnc(three_period_cost_two_classes, [0.5, 0.5, 0.5],
                                  args=(x0, u_bar, self.h_cost, self.r_time[sim], mu, norm_lbda),
                                  messages=0, approx_grad=True, bounds=[(0, 1), (0, 1), (0, 1)])
            u1  = np.around(new_result[0][0]*self.N)
            return [u1, self.N-u1]

    def _get_new_alloc_three_period_d2_tnc(self, sim, safe_coef):
        mu = [1 / x for x in self.w_mu]
        safety = [safe_coef * self.N ** .5, 0]
        x0 = [max(
            (self.ward_alloc[sim][i] - self.ward_capac[sim][i] + self.queue_length[sim][i] - safety[i]) / float(self.N),
            0)
              for i in range(self.k)]
        norm_lbda = [1 / (x * float(self.N)) for x in self.l_arr]
        u_bar = [(x0[i] + norm_lbda[i] * self.r_time[sim]) / (1 + mu[i] * self.r_time[sim]) for i in range(self.k)]
        if x0[0] + x0[1] < 1:
            return self.dedicated_alloc[sim]
        elif x0[0] >= (mu[0] - norm_lbda[0]) * self.r_time[0] + 1:
            return [self.N, 0]
        else:
            new_result = fmin_tnc(three_period_cost_two_classes, [0.5, 0.5, 0.5],
                                  args=(x0, u_bar, self.h_cost, self.r_time[sim], mu, norm_lbda),
                                  messages=0, approx_grad=True, bounds=[(0, 1), (0, 1), (0, 1)])
            u1 = np.around(new_result[0][0] * self.N)
            return [u1, self.N - u1]

    def _get_new_alloc_three_period_d3_tnc(self, sim, safe_coef):
        mu = [1 / x for x in self.w_mu]
        safety = [safe_coef * self.N ** .5, 0]
        x0 = [max(
            (self.ward_alloc[sim][i] - self.ward_capac[sim][i] + self.queue_length[sim][i] - safety[i]) / float(self.N),
            0)
              for i in range(self.k)]
        norm_lbda = [1 / (x * float(self.N)) for x in self.l_arr]
        rho = [norm_lbda[i] / mu[i] for i in range(self.k)]
        u_bar = [(x0[i] + norm_lbda[i] * self.r_time[sim]) / (1 + mu[i] * self.r_time[sim]) for i in range(self.k)]
        if x0[0] <= rho[0] and x0[1] <= rho[1]:
            return self.dedicated_alloc[sim]
        elif x0[0] >= (mu[0] - norm_lbda[0]) * self.r_time[0] + 1:
            return [self.N, 0]
        else:
            new_result = fmin_tnc(three_period_cost_two_classes, [0.5, 0.5, 0.5],
                                  args=(x0, u_bar, self.h_cost, self.r_time[sim], mu, norm_lbda),
                                  messages=0, approx_grad=True, bounds=[(0, 1), (0, 1), (0, 1)])
            u1 = np.around(new_result[0][0] * self.N)
            if sum(x0) < 1:
                if u1 < np.floor(rho[0]*self.N) + 1:
                    return [np.floor(rho[0]*self.N) + 1, self.N - (np.floor(rho[0]*self.N) + 1)]
                elif self.N - u1 < np.floor(rho[1]*self.N) + 1:
                    return [self.N - (np.floor(rho[1] * self.N) + 1), np.floor(rho[1] * self.N) + 1]
                else:
                    return [u1, self.N - u1]
            else:
                return [u1, self.N - u1]

    def _get_new_alloc_three_period_d2_slsqp(self, sim, safe_coef, rnd):
        mu = [1 / x for x in self.w_mu]
        safety = [safe_coef * self.N ** .5, 0]
        x0 = [max(
            (self.ward_alloc[sim][i] - self.ward_capac[sim][i] + self.queue_length[sim][i] - safety[i]) / float(self.N),
            0)
              for i in range(self.k)]
        norm_lbda = [1 / (x * float(self.N)) for x in self.l_arr]
        u_bar = [(x0[i] + norm_lbda[i] * self.r_time[sim]) / (1 + mu[i] * self.r_time[sim]) for i in range(self.k)]
        if x0[0] + x0[1] <= 1:
            return self.dedicated_alloc[sim]
        elif x0[0] >= (mu[0] - norm_lbda[0]) * self.r_time[0] + 1:
            return [self.N, 0]
        else:
            cons = ({'type': 'ineq', 'fun': lambda x: -x[0] + 1},
                    {'type': 'ineq', 'fun': lambda x: -x[1] + 1},
                    {'type': 'ineq', 'fun': lambda x: -x[2] + 1},
                    {'type': 'ineq', 'fun': lambda x: -x[3] + 1},
                    {'type': 'ineq', 'fun': lambda x: -x[4] + 1},
                    {'type': 'ineq', 'fun': lambda x: -x[5] + 1},
                    {'type': 'ineq', 'fun': lambda x: -(x[0] + x[3]) + 1},
                    {'type': 'ineq', 'fun': lambda x: -(x[1] + x[4]) + 1},
                    {'type': 'ineq', 'fun': lambda x: -(x[2] + x[5]) + 1}
                    )
            new_result = minimize(three_period_cost_two_classes_0, [0.5, 0.5, 0.5, 0.5, 0.5, 0.5],
                              args=(x0, u_bar, self.h_cost, self.r_time[sim], mu, norm_lbda), jac=False,
                              constraints=cons, method='SLSQP', options={'disp': False, 'ftol': 1e-08})
            if rnd == 0:
                new_alloc = [np.around(new_result['x'][0] * self.N), np.around(new_result['x'][3] * self.N)]
                new_alloc = rounding_pref_2_class(new_alloc, self.N)
                return new_alloc
            elif rnd == 1:
                u1 = np.around(new_result['x'][0] * self.N)
                return [u1, self.N - u1]
            else:
                print 'error'

    def _get_new_alloc_three_period_d1_slsqp(self, sim, safe_coef, rnd):
        mu = [1 / x for x in self.w_mu]
        safety = [safe_coef * self.N ** .5, 0]
        x0 = [max(
            (self.ward_alloc[sim][i] - self.ward_capac[sim][i] + self.queue_length[sim][i] - safety[i]) / float(self.N),
            0)
              for i in range(self.k)]
        norm_lbda = [1 / (x * float(self.N)) for x in self.l_arr]
        rho = [norm_lbda[i] / mu[i] for i in range(self.k)]
        u_bar = [(x0[i] + norm_lbda[i] * self.r_time[sim]) / (1 + mu[i] * self.r_time[sim]) for i in range(self.k)]
        if x0[0] <= rho[0] and x0[1] <= rho[1]:
            return self.dedicated_alloc[sim]
        elif x0[0] >= (mu[0] - norm_lbda[0]) * self.r_time[0] + 1:
            return [self.N, 0]
        else:
            cons = ({'type': 'ineq', 'fun': lambda x: -x[0] + 1},
                    {'type': 'ineq', 'fun': lambda x: -x[1] + 1},
                    {'type': 'ineq', 'fun': lambda x: -x[2] + 1},
                    {'type': 'ineq', 'fun': lambda x: -x[3] + 1},
                    {'type': 'ineq', 'fun': lambda x: -x[4] + 1},
                    {'type': 'ineq', 'fun': lambda x: -x[5] + 1},
                    {'type': 'ineq', 'fun': lambda x: -(x[0] + x[3]) + 1},
                    {'type': 'ineq', 'fun': lambda x: -(x[1] + x[4]) + 1},
                    {'type': 'ineq', 'fun': lambda x: -(x[2] + x[5]) + 1}
                    )
            new_result = minimize(three_period_cost_two_classes_0, [0.5, 0.5, 0.5, 0.5, 0.5, 0.5],
                                  args=(x0, u_bar, self.h_cost, self.r_time[sim], mu, norm_lbda), jac=False,
                                  constraints=cons, method='SLSQP', options={'disp': False, 'ftol': 1e-08})
            if rnd == 0:
                new_alloc = [np.around(new_result['x'][0] * self.N), np.around(new_result['x'][3] * self.N)]
                new_alloc = rounding_pref_2_class(new_alloc, self.N)
                return new_alloc
            elif rnd == 1:
                u1 = np.around(new_result['x'][0] * self.N)
                return [u1, self.N - u1]
            else:
                print 'error'

    def _get_new_alloc_three_period(self, sim):
        mu = [1 / x for x in self.w_mu]
        # safety = [self.N ** .5, 0]
        safety = [0,0]
        x0 = [max(
            (self.ward_alloc[sim][i] - self.ward_capac[sim][i] + self.queue_length[sim][i] - safety[i]) / float(
                self.N), 0)
              for i in range(self.k)]
        norm_lbda = [1 / (x * float(self.N)) for x in self.l_arr]
        rho = [norm_lbda[i] / mu[i] for i in range(self.k)]
        u_bar = [(x0[i] + norm_lbda[i] * self.r_time[sim]) / (1 + mu[i] * self.r_time[sim]) for i in range(self.k)]
        if x0[0] <= rho[0] and x0[1] <= rho[1]:
            # if x0[0] + x0[1] < 1:
            return self.dedicated_alloc[sim]
        elif x0[0] >= (mu[0] - norm_lbda[0]) * self.r_time[0] + 1:
            return [self.N, 0]
        else:
            cons = ({'type': 'ineq', 'fun': lambda x: -x[0] + 1},
                    {'type': 'ineq', 'fun': lambda x: -x[1] + 1},
                    {'type': 'ineq', 'fun': lambda x: -x[2] + 1},
                    {'type': 'ineq', 'fun': lambda x: -x[3] + 1},
                    {'type': 'ineq', 'fun': lambda x: -x[4] + 1},
                    {'type': 'ineq', 'fun': lambda x: -x[5] + 1},
                    {'type': 'ineq', 'fun': lambda x: -(x[0] + x[3]) + 1},
                    {'type': 'ineq', 'fun': lambda x: -(x[1] + x[4]) + 1},
                    {'type': 'ineq', 'fun': lambda x: -(x[2] + x[5]) + 1}
                    )
            new_result = minimize(three_period_cost_two_classes_0, [0.5, 0.5, 0.5, 0.5, 0.5, 0.5],
                                  args = (x0, u_bar, self.h_cost, self.r_time[sim], mu, norm_lbda), jac = False,
                                  constraints = cons, method = 'SLSQP', options = {'disp': False,'ftol': 1e-08})
            # new_alloc = [np.around(new_result['x'][0] * self.N), np.around(new_result['x'][3] * self.N)]
            # new_alloc = rounding_pref_2_class(new_alloc, self.N)
            u1 = np.around(new_result['x'][0] * self.N)
            return [u1, self.N - u1]

# ------------------- Numerical Methods -------------------
def rounding_pref_2_class(solution, N):
    new_solution = solution
    while sum(new_solution) != N:
        if sum(new_solution) > N:
            if new_solution[1] > 0:
                new_solution[1] -= 1
            else:
                new_solution[0] -= 1
        else:
            new_solution[0] += 1

    return new_solution

def mean_confidence_interval(data, confidence=0.95):
    a = 1.0*np.array(data)
    n = len(a)
    m, se = np.mean(a), scipy.stats.sem(a)
    h = se * sp.stats.t._ppf((1+confidence)/2., n-1)
    return m, m-h, m+h

def vs_exp(l, period, a, t, vary):
    if vary:
        original_t = t
        curr_t = t
        lv = 2*a*float(l)
        U = 1.1
        while(U > lv/(2*a*float(l))):
            U = np.random.uniform()
            curr_t = curr_t + np.random.exponential(1/float(l*a*2))
            lv = l*a + a*l*np.sin(2*np.pi*curr_t/period)
        return (curr_t - original_t)
    else:
        return np.random.exponential(1/l)

def ode_sys(x, t, y_0, t_0, lbda, mu, cnt):
    tot = []
    for i in range(cnt):
        tot.append(solve_ode_sys(t, y_0, t_0, lbda, mu, i))
    return lbda[cnt]-mu[cnt]*min(x, max(0,1-sum(tot)))

def solve_ode_sys(t, y_0, t_0, lbda, mu, cnt):
    y0 = y_0[cnt]
    t0 = t_0[cnt]
    tF = t
    tt = np.linspace(t0, tF, 2)
    yy = odeint(ode_sys, y0, tt, args=(y_0, t_0, lbda, mu, cnt))
    return yy[1][0]

def ode_sys_complete(t, y0, t0, lbda, mu, cnt):
    alloc = []
    for i in range(cnt):
        sol = solve_ode_sys(t, y0, t0, lbda, mu, i)
        int_sol = (sol - y0[i] - lbda[i]*t)/-mu[i]
        alloc.append(int_sol)
    return alloc

def two_class_cost_foc(u, x, mu_bar, cost, t, serv, arr):
    rho = arr / float(serv)
    if x >= rho:
        if u <= mu_bar:
            return -cost*t-cost*serv*(t**2)/2.0
        elif mu_bar < u < x:
            sig = (x-u)/(u*serv-arr)
            return -cost*sig-cost*serv*(sig**2)/2.0
        elif u >= x:
            return 0
    else:
        if u < x:
            return -cost*t-cost*serv*t**2/2.0
        elif x <= u < rho:
            nu = (-1/serv)*np.log((rho-u)/(rho-x))
            return -cost*max(t-nu,0)-cost*serv*(max(0,t-nu)**2)/2.0
        elif u >= rho:
            return 0

def two_class_cost_foc_combined(u, x, u_bar, cost, tau, mu, lbda):
    a = two_class_cost_foc(u, x[0], u_bar[0], cost[0], tau, mu[0], lbda[0])
    b = two_class_cost_foc(1-u, x[1],u_bar[1],cost[1], tau, mu[1], lbda[1])
    return a - b

def find_foc_root(x_bar, u_bar, cost, tau, mu, lbda):
    x0 = two_class_cost_foc_combined(0, x_bar, u_bar, cost, tau, mu, lbda)
    x1 = two_class_cost_foc_combined(1, x_bar, u_bar, cost, tau, mu, lbda)
    if x0 < 0 and x1 < 0:
        return 1
    elif x0 > 0 and x1 > 0:
        return 0
    else:
        return brentq(two_class_cost_foc_combined, 0, 1, args=(x_bar, u_bar, cost, tau, mu, lbda))

def two_period_cost(u_1, u_2, x, u_bar, cost, tau, mu, lbda):
    x_2 = get_trajectory(u_1, x, mu, lbda, tau, tau)
    u_bar_2 = (x_2+lbda*tau)/(1+mu*tau)
    return class_cost(u_1, x, u_bar, cost, tau, mu, lbda) + \
           class_cost(u_2, x_2, u_bar_2, cost, tau, mu, lbda)

def two_period_cost_two_classes(u, x_bar, u_bar, cost, tau, mu, lbda_bar):
    return two_period_cost(u[0], u[1], x_bar[0], u_bar[0], cost[0], tau, mu[0], lbda_bar[0]) +\
           two_period_cost(1-u[0],  1-u[1], x_bar[1], u_bar[1], cost[1], tau, mu[1], lbda_bar[1])

def three_period_cost(u_1, u_2, u_3, x, u_bar, cost, tau, mu, lbda):
    x_2 = get_trajectory(u_1, x, mu, lbda, tau, tau)
    x_3 = get_trajectory(u_2, x_2, mu, lbda, tau, tau)
    u_bar_2 = (x_2+lbda*tau)/(1+mu*tau)
    u_bar_3 = (x_3 + lbda * tau) / (1 + mu * tau)
    return class_cost(u_1, x, u_bar, cost, tau, mu, lbda) + class_cost(u_2, x_2, u_bar_2, cost, tau, mu, lbda) +\
           class_cost(u_3, x_3, u_bar_3, cost, tau, mu, lbda)

def three_period_cost_two_classes_0(u, x_bar, u_bar, cost, tau, mu, lbda_bar):
    return three_period_cost(u[0], u[1], u[2], x_bar[0], u_bar[0], cost[0], tau, mu[0], lbda_bar[0]) + \
           three_period_cost(u[3], u[4], u[5], x_bar[1], u_bar[1], cost[1], tau, mu[1], lbda_bar[1])

def three_period_cost_two_classes(u, x_bar, u_bar, cost, tau, mu, lbda_bar):
    return three_period_cost(u[0], u[1], u[2], x_bar[0], u_bar[0], cost[0], tau, mu[0], lbda_bar[0]) + \
           three_period_cost(1-u[0], 1-u[1], 1-u[2], x_bar[1], u_bar[1], cost[1], tau, mu[1], lbda_bar[1])

def m_period_cost(u_arr, m, x, u_bar, cost, tau, mu, lbda):
    total_cost = 0
    traj = x
    for i in range(m):
        total_cost += class_cost(u_arr[i], traj, u_bar, cost, tau, mu, lbda)
        traj = get_trajectory(u_arr[i], traj, mu, lbda, tau, tau)
    return total_cost

def m_period_cost_two_classes(u, m, x_bar, u_bar, cost, tau, mu, lbda_bar):
    total_cost = 0
    for i in range(2):
        total_cost += m_period_cost(u[m*i:m*i+m])
    return total_cost

def class_cost(u, x, u_bar, cost, tau, mu, lbda):
    rho = lbda/float(mu)
    if x > rho:
        if u_bar >= u:
            return cost * (x - u) * tau + cost * (lbda - mu * u) * (tau ** 2) / 2.0
        elif u_bar < u < x:
            return (cost * (x - u) ** 2) / (2 * (mu * u - lbda))
        elif u >= x:
            return 0
    else:
        if x > u:
            return cost * (x - u) * tau + cost * (lbda - mu * u) * (tau ** 2) / 2.0
        elif x <= u < rho:
            nu = -np.log((rho - u)/(rho - x))/mu
            return (cost * (lbda - u * mu) * (max(0,tau - nu)) ** 2) / 2
        elif u >= rho:
            return 0

def get_trajectory(u, x, mu, lbda, tau, t):
    rho = lbda / float(mu)
    if x >= u:
        if x >= u + (u * mu - lbda) * tau:
            return x + (lbda - u*mu)*t
        else:
            sig = (x - u)/(u*mu - lbda)
            if t > sig:
                return rho + np.exp(-mu*(t-sig))*(u-rho)
            else:
                return x + (lbda - u * mu)*t
    else:
        if u >= rho:
            return rho + np.exp(-mu*t)*(x - rho)
        else:
            nu = -np.log((rho - u)/(rho - x))/mu
            if nu >= tau:
                return rho + np.exp(-mu*t)*(x - rho)
            else:
                return u + (lbda - u*mu)* (t - nu)

def ErlangC(p, n):
    _sum = 0
    for i in range(n):
        _sum += (n*p)**i/factorial(i)
    return  (n*p)**n/factorial(n)/(_sum*(1-p) + (n*p)**n/factorial(n))

def queue_length(lbda, mu, n, N):
    p = lbda/(mu*n)
    return ErlangC(p, n)*(p)/(1-p)

def cost_func(u, in_lbda, in_mu, N, cost):
    _sum = 0
    lbda = [1 / float(x) for x in in_lbda]
    mu = [1 / float(x) for x in in_mu]
    for ind, c in enumerate(cost):
        _sum += c*queue_length(lbda[ind], mu[ind], u[ind], N)
    return _sum

def find_dedicated_2_class(N, in_lbda, in_mu, cost):
    bound = []
    lbda = [1/float(x) for x in in_lbda]
    mu = [1/float(x) for x in in_mu]
    for i in range(len(lbda)):
        bound.append(lbda[i]/mu[i])
    min_val = float('inf')
    min_u = 0
    for i in range(N+1):
        if bound[0] < i and bound[1] < N - i:
            u = [i, N-i]
            cost_val = cost_func(u, in_lbda, in_mu, N, cost)
            if cost_val < min_val:
                min_u = i
                min_val = cost_val
    return [min_u, N-min_u]

def writeLog(fil, table):
    c1 = csv.writer(fil)
    for val in table:
        c1.writerow(val)

# Modify simulation below:
if __name__ == "__main__":
    # ================ Input Variables ================
    Total_Time = 72
    # Scale by nurses
    Nurses = 50
    # lbda_out = [1.0/(.24*Nurses), 1.0/(.24*Nurses)]
    lbda_out = [1/23.0, 1/23.0]
    # mu_out = [1.0/.5, 1.0/.5]
    mu_out = [1, 1]
    std_out = [1, 1]
    theta_out = [10000, 10000]
    tau_out = [24, 24, 24]
    k_out = 2
    hcost_out = [2,1]
    q_cap_out = [float('inf'), float('inf')]
    # Parallel simulation variables
    tot_par = 3
    s_alloc_out = [[Nurses/2,Nurses/2], [Nurses,Nurses], [Nurses/2, Nurses/2]]
    rebalance1 = [4, 1, 0]
    cont_out = [0, 0, 1]
    preemption_out = [0, 0, 1]
    time_vary = False
    # Trial variables
    trials = 500

    print 'Nurses: ' + str(Nurses)
    print 'tau: ' + str(tau_out)
    dataset_arr, dataset_hc, dataset_st = [[] for x in range(tot_par)], [[] for x in range(tot_par)], [[] for x in range(tot_par)]
    dataset_wq, dataset_ww = [[[] for y in range(k_out)] for x in range(tot_par)], [[[] for y in range(k_out)] for x in range(tot_par)]
    for t in range(trials):
        start_time = time.clock()
        s = Simulation(Total_Time, Nurses, lbda_out, mu_out, std_out, theta_out, tau_out, k_out, hcost_out, q_cap_out,
                       s_alloc_out, tot_par, rebalance1, cont_out, preemption_out, time_vary)
        # class0 = [0 for x in range(50)]
        # class1 = [1 for x in range(50)]
        # class0mu = [mu_out[0] for x in range(50)]
        # class1mu = [mu_out[1] for x in range(50)]
        # class0ab = [theta_out[0] for x in range(50)]
        # class1ab = [theta_out[1] for x in range(50)]
        # s.set_preexisting_rnd(class0+class1,class0mu+class1mu,class0ab+class1ab)
        # Choose Option to utilize time-varying arrivals here
        s.generate_arrivals(time_vary)
        s.simulate(False, True)
        # Save data
        for p in range(tot_par):
            dataset_arr[p].append(s.arrival_count[p]/float(s.t[p]))
            dataset_hc[p].append(s.holding_cost[p]/float(s.t[p]))
            dataset_st[p].append(s.time_server_occupied[p]/(Nurses*float(s.t[p])))
            for c in range(k_out):
                dataset_wq[p][c].append(s.weighted_queue[p][c]/float(s.t[p]))
                dataset_ww[p][c].append(s.weighted_ward[p][c]/float(s.t[p]))
        print 'finished: ' + str(t)
        print str(time.clock() - start_time) + ' secs'
    for p in range(tot_par):
        print "\nSimulation " + str(p)
        print "Arrivals CI: " + str(mean_confidence_interval(dataset_arr[p]))
        print "Holding Cost CI: " + str(mean_confidence_interval(dataset_hc[p]))
        # print "Holding Cost" + str(dataset_hc[p])
        print "Server Time CI: " + str(mean_confidence_interval(dataset_st[p]))
        for i in range(k_out):
            print "Queue length for ward CI " + str(i) + ": " + str(mean_confidence_interval(dataset_wq[p][i]))
            # print "Queue length" + str(dataset_wq[p][i])
            print "Headcount for ward CI" + str(i) + ": " + str(mean_confidence_interval(dataset_ww[p][i]))
            # print "Headcount for ward CI" + str(dataset_ww[p][i])
        print ' '

    #
    #
    # fil0 = open(os.getcwd() + "/Sim_Rebalance_6.csv", "wb")
    # fil1 = open(os.getcwd() + "/Sim_No_Rebalance_6.csv", "wb")
    #
    # writeLog(fil0, s.statistics[0])
    # writeLog(fil1, s.statistics[1])
    #
