import numpy as np 
import heapq
import os
import time
import random 
import csv
import scipy as sp
import scipy.stats

# Global Variables for easier use in the simulation.




# ----------------------------------- Parameters -----------------------------------
pm = 0				# Number of parallel simulations
k = 0 				# Number of patient types
l_arr = [] 			# Array containig Arrival rate for each simulation
l_aban = []			# Array containing abandonment rate for each simulation
w_mu = [] 			# Array containing mean service time
w_std = [] 			# Array containing STD of service time
h_cost = []			# Array containing holding cost of each ward queue
n_free = [] 		# Array containing number of nurses free
preempt = []		# Array determining whether or not simulation is running preemption

# ----------------------------------- Simulation Variables -----------------------------------
t = [] 				# Array to hold current time for each simulation
r_time_arr = [] 	# Array holding next time of rebalance for each simulation
queue = [] 			# Array containing the patients in the queues of each ward for each simulation
q_capac = [] 		# Array containing the capacity for each queue
Events = [] 		# List containing Events where patients leave the ward
A_Events = [] 		# List containing Events of arrival for each patient type

a_queue_count = [] 			# Counter for number of preloaded arrivals

queue_length = [] 			# tracks queue length for each ward
ward_alloc = [] 			# tracks current allocation of nurses in each ward
ward_capac = [] 			# tracks current capacity available for each ward
ward_nurse_def = []			# tracks number of nurses ward still needs to receive after rebalance
ward_assignment = [] 		# tracks wards that still need to receive nurses after rebalance

last_arrival = []			# tracks last patient that arrived in a ward

# ---------------------- Counters for the simulation ----------------------
Arrival_Count = [] 			# Counts total number of arrivals
Abandonment_Count = [] 		# Counts total number of abandoments
Balk_Count = [] 			# Counter for patients who balk
Treated = [] 				# Counter for number of patients treated
holding_cost = [] 			# Counter for holding cost
time_server_occupied = [] 	# Counter for time servers are occupied
weighted_queue = [] 		# Counter to calculate average queue length
weighted_ward = []			# Counter to calculate average headcount in ward

Counter_Record = []

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

def mean_confidence_interval(data, confidence=0.95):
    
    a = 1.0*np.array(data)
    n = len(a)
    m, se = np.mean(a), scipy.stats.sem(a)
    h = se * sp.stats.t._ppf((1+confidence)/2., n-1)
    return m, m-h, m+h

# Used to reset simulation parameters for new trial
def resetvar():

	global pm, k, l_arr, l_aban, w_mu, w_std, h_cost, n_free, preempt, t, r_time_arr, queue, q_capac, Events
	global a_queue_count, queue_length, ward_alloc, ward_capac, ward_nurse_def, ward_assignment, Arrival_Count
	global Balk_Count, Treated, holding_cost, time_server_occupied, Abandonment_Count, A_Events, weighted_ward, weighted_queue
	global Counter_Record


	new_row = []
	new_row.append(Arrival_Count)
	new_row.append(holding_cost)
	new_row.append(time_server_occupied)
	new_row.append(weighted_queue)
	new_row.append(weighted_ward)
	new_row.append(t)

	Counter_Record.append(new_row)

	# ----------------------------------- Parameters -----------------------------------
	pm = 0				# Number of parallel simulations
	k = 0 				# Number of patient types
	l_arr = [] 			# Array containig Arrival rate for each simulation
	l_aban = []			# Array containing abandonment rate for each simulation
	w_mu = [] 			# Array containing mean service time
	w_std = [] 			# Array containing STD of service time
	h_cost = []			# Array containing holding cost of each ward queue
	n_free = [] 		# Array containing number of nurses free
	preempt = []		# Array determining whether or not simulation is running preemption

	# ----------------------------------- Simulation Variables -----------------------------------
	t = [] 				# Array to hold current time for each simulation
	r_time_arr = [] 	# Array holding next time of rebalance for each simulation
	queue = [] 			# Array containing the patients in the queues of each ward for each simulation
	q_capac = [] 		# Array containing the capacity for each queue
	Events = [] 		# List containing Events where patients leave the ward
	A_Events = [] 		# List containing Events of arrival for each patient type

	a_queue_count = [] 			# Counter for number of preloaded arrivals

	queue_length = [] 			# tracks queue length for each ward
	ward_alloc = [] 			# tracks current allocation of nurses in each ward
	ward_capac = [] 			# tracks current capacity available for each ward
	ward_nurse_def = []			# tracks number of nurses ward still needs to receive after rebalance
	ward_assignment = [] 		# tracks wards that still need to receive nurses after rebalance


	# ---------------------- Counters for the simulation ----------------------
	Arrival_Count = [] 			# Counts total number of arrivals
	Abandonment_Count = [] 		# Counts total number of abandoments
	Balk_Count = [] 			# Counter for patients who balk
	Treated = [] 				# Counter for number of patients treated
	holding_cost = [] 			# Counter for holding cost
	time_server_occupied = [] 	# Counter for time servers are occupied
	weighted_queue = [] 		# Counter to calculate average queue length
	weighted_ward = []			# Counter to calculate average headcount in ward


def rebalance(N, sim):

	global queue_length

	global ward_capac
	global ward_alloc
	global ward_assignment
	global ward_nurse_def
	global Events
	global A_Events
	global n_free 
	global pm
	global k
	global t
	global preempt

	# print queue_length[sim]
	# print ward_capac[sim]
	# print n_free

	# Set variable to hold old allocation
	old_alloc = [0]*k


	# -------------------- START: Procedure for setting new allocation --------------------
	total = max(sum(queue_length[sim]), k)
	s = 0
	for i in range(k-1):
		old_alloc[i] = ward_alloc[sim][i]
		ward_alloc[sim][i] = queue_length[sim][i]*(N-k)/total + 1
		ward_capac[sim][i] = ward_capac[sim][i] + ward_alloc[sim][i] - old_alloc[i]
		s += queue_length[sim][i]*(N-k)/total

	# -------------------- END: Procedure for setting new allocation --------------------

	# Setting new allocation and calculating new capacities
	old_alloc[k-1] = ward_alloc[sim][k-1]
	ward_alloc[sim][k-1] = (N - k) - s + 1
	ward_capac[sim][k-1] = ward_capac[sim][k-1] + ward_alloc[sim][k-1] - old_alloc[k-1]
	
	# Procedure to output information to check for consistency
	# print 'phase i'
	# print ward_capac[sim]
	# print ward_alloc[sim]
	# print n_free[sim]
	
	# Preemption procedure
	requeue = []
	if preempt[sim] == 1:
	
		for i in range(k):
			requeue.append([])

		# Changes the ward patients back to queue patients after readjusting service time
		# to remaining service time. 
		for pat in Events[sim]:
			remaining_service = pat.get_time() - t[sim]
			pat.set_serv(remaining_service)
			pat.set_location('abandonment')
			requeue[pat.get_pt()].append(pat)
		
		# Reset the wards so they are empty and nurses are all freed
		Events[sim] = []
		ward_capac[sim] = [x for x in ward_alloc[sim]]
		n_free[sim] = N

		# Procecure pushes all patients in wards back to their respective queues
		# We then redistribute the patients later in the next block of code
		for i in range(k):
			incoming_length = len(requeue[i])
			queue[sim][i] = requeue[i] + queue[sim][i]
			queue_length[sim][i] += incoming_length



	# For each type of class of patients we will check to see if we can refill the wards
	for i in range(k):
		# If there are nurses free and the ward has capacity and has a queue, we assign patients
		while(n_free[sim] > 0 and ward_capac[sim][i] > 0 and queue_length[sim][i] > 0):
			# Get next patient in the queue
			next_patient = queue[sim][i].pop(0)
			queue_length[sim][i] -= 1
			
			# Change patient type to ward patient from abandoner
			next_patient.set_time(t[sim] + next_patient.get_serv())
			next_patient.set_location('ward')
			
			# Add patient to Events list
			heapq.heappush(Events[sim], next_patient)

			# Update counters 
			n_free[sim] -= 1
			ward_capac[sim][i] -= 1

	# Determines deficits in wards by checking for excess capacity when no nurses are free		
	for i in range(k):
		if n_free[sim] == 0 and ward_capac[sim][i] > 0:
			ward_nurse_def[sim][i] = -ward_capac[sim][i]
			ward_assignment[sim].append(i)
	
	# Procedure to output information to check for consistency
	# print 'phase ii'
	# print ward_capac[sim]
	# print ward_alloc[sim]
	# print n_free[sim]

# --------------------------------- Arrival Event for continuous policy ---------------------------------
def arrival_event_cont(event, sim):
	global Arrival_Count
	global Balk_Count
	global ward_capac
	global w_mu
	global w_std
	global t
	global l_aban
	global Events
	global A_Events
	global l_arr
	global queue 
	global n_free
	global q_capac
	global a_queue_count
	global pm
	global queue_length
	global preempt
	global last_arrival
	global h_cost

	Arrival_Count[sim] += 1

	pt = event.get_pt()

	if n_free[sim] > 0:
		# Get service and abandonment times for the patient
		new_serv_time = event.get_serv()
		new_aban = event.get_aban()

		# Create new patient with updated event time for ward
		new_ward_arr = Patient(new_serv_time + t[sim], t[sim], pt, 'ward', 
			new_aban, new_serv_time)

		# Push new patient to Event
		heapq.heappush(Events[sim], new_ward_arr)

		# Adjust the number of nurses free
		ward_capac[sim][pt] = ward_capac[sim][pt] - 1
		n_free[sim] = n_free[sim] - 1

		# Set last ward arrival 
		# last_arrival[sim][pt] = new_ward_arr

	else:
		min_cost = h_cost[pt]
		min_pt = pt
		

		if preempt[sim] == 1:
			
			for h in range(k):
				if h_cost[h] < min_cost and (ward_alloc[sim][h]-ward_capac[sim][h]) > 0:
					min_cost = h_cost[h]
					min_pt = h

		if h_cost[pt] > min_cost:
			for pat_ind, pat_var in reversed(list(enumerate(Events[sim]))):
				
				if pat_var.get_pt() == min_pt:
					pat1 = pat_var
					del Events[sim][pat_ind]
					heapq.heapify(Events[sim])
					
					# print n_free[sim]
					break


			# Remove last patient that arrived from the ward with lowest cost
			# pat1 = last_arrival[sim][min_pt]
			
			pat1.set_serv(pat1.get_time()-t[sim])
			pat1.set_location('abandonment')
			queue[sim][min_pt].insert(0,pat1)
			queue_length[sim][min_pt] += 1
			heapq.heapify(Events[sim])

			# Get service and abandonment times for the patient
			new_serv_time = event.get_serv()
			new_aban = event.get_aban()

			# Create new patient with updated event time for ward
			new_ward_arr = Patient(new_serv_time + t[sim], t[sim], pt, 'ward', 
				new_aban, new_serv_time)

			# Push new patient to Event
			heapq.heappush(Events[sim], new_ward_arr)

			# Adjust the number of nurses free
			ward_capac[sim][min_pt] = ward_capac[sim][min_pt] + 1
			ward_capac[sim][pt] = ward_capac[sim][pt] - 1


			# Set last ward arrival
			# last_arrival[sim][min_pt] = new_ward_arr

		elif len(queue[sim][pt]) < q_capac[pt]:
			# Get service and abandonment times for the patient
			new_arrival_time = event.get_serv()
			new_aban = event.get_aban()
			
			# Create new patient with updated event time for queue
			new_queue_arr = Patient(new_aban + t[sim], t[sim], pt, 'abandonment', 
				new_aban, new_arrival_time)
			
			# Push new patient to Event
			queue[sim][pt].append(new_queue_arr)
			queue_length[sim][pt] += 1

		#Case where queue is full and patient leaves system
		else:
			Balk_Count[sim] += 1

	# Remove Patient Arrival Event and reduce count
	if t[sim] == 0:
		heapq.heappop(Events[sim])
	else: 
		heapq.heappop(A_Events[sim])
		a_queue_count[sim][pt] -= 1

	#Create new arrival and replace old arrival
	if a_queue_count[sim][pt] == 0 and t[sim] != 0:

		# Generate new arrival time for patient and create new patient
		new_arrival_time = np.random.exponential(l_arr[pt])
		new_arr = Patient(new_arrival_time + t[sim], t[sim], pt, 'arrival', 
			np.random.exponential(l_aban[pt]), np.random.exponential(w_mu[pt]))
		# Add new arrival to all queues
		for i in range(pm):

			heapq.heappush(A_Events[i], new_arr)
			a_queue_count[i][pt] += 1


# --------------------------------- Departure Event for continuous policy ---------------------------------
def departure_event_cont(event, sim):
	global Events
	global A_Events
	global queue
	global w_mu
	global w_std
	global ward_capac
	global t
	global n_free
	global Treated
	global queue_length
	global ward_assignment
	global ward_nurse_def
	global k 
	global h_cost
	global last_arrival

	pt = event.get_pt()


	max_cost = 0
	max_pt = pt
	for h in range(k):
		if h_cost[h] >= max_cost and queue_length[sim][h] > 0:
			max_cost = h_cost[h]
			max_pt = h

	if queue_length[sim][max_pt] > 0:
		# Get next patient and remove from Event list (abandoner) to change to a ward patient
		next_patient = queue[sim][max_pt].pop(0)
		# Events[sim].remove(next_patient)
		queue_length[sim][max_pt] -= 1
		ward_capac[sim][pt] = ward_capac[sim][pt] + 1
		ward_capac[sim][max_pt] = ward_capac[sim][max_pt] - 1

		
		# Change patient type to ward patient from abandoner
		next_patient.set_time(t[sim] + next_patient.get_serv())
		next_patient.set_location('ward')
		
		# Heapify the event list and replace the last departure with the new patient
		heapq.heappop(Events[sim])
		heapq.heappush(Events[sim], next_patient)
	else:
		heapq.heappop(Events[sim])
		n_free[sim] += 1
		ward_capac[sim][pt] += 1

	# if event == last_arrival[sim][pt]:
	# 	print 'check'
	# 	print event

	Treated[sim] += 1

# --------------------------------- Arrival Event for non-continuous policies ---------------------------------
def arrival_event(event, sim):
	global Arrival_Count
	global Balk_Count
	global ward_capac
	global w_mu
	global w_std
	global t
	global l_aban
	global Events
	global A_Events
	global l_arr
	global queue 
	global n_free
	global q_capac
	global a_queue_count
	global pm
	global queue_length

	Arrival_Count[sim] += 1

	pt = event.get_pt()
	
	#Case where designated ward is open
	if ward_capac[sim][pt] > 0 and n_free[sim] > 0:
		# Get service and abandonment times for the patient
		new_serv_time = event.get_serv()
		new_aban = event.get_aban()

		# Create new patient with updated event time for ward
		new_ward_arr = Patient(new_serv_time + t[sim], t[sim], pt, 'ward', 
			new_aban, new_serv_time)

		# Push new patient to Event
		heapq.heappush(Events[sim], new_ward_arr)

		# Adjust the capacities of the wards and nurses free
		ward_capac[sim][pt] = ward_capac[sim][pt] - 1
		n_free[sim] = n_free[sim] - 1
		
	#Case where ward is full and patient is sent into queue
	elif len(queue[sim][pt]) < q_capac[pt]:
		# Get service and abandonment times for the patient
		new_arrival_time = event.get_serv()
		new_aban = event.get_aban()
		
		# Create new patient with updated event time for queue
		new_queue_arr = Patient(new_aban + t[sim], t[sim], pt, 'abandonment', 
			new_aban, new_arrival_time)
		
		# Push new patient to Event
		queue[sim][pt].append(new_queue_arr)
		queue_length[sim][pt] += 1

	#Case where queue is full and patient leaves system
	else:
		Balk_Count[sim] += 1

	# Remove Patient Arrival Event and reduce count
	if t[sim] == 0:
		heapq.heappop(Events[sim])
	else: 
		heapq.heappop(A_Events[sim])
		a_queue_count[sim][pt] -= 1

	#Create new arrival and replace old arrival
	if a_queue_count[sim][pt] == 0 and t[sim] != 0:

		# Generate new arrival time for patient and create new patient
		new_arrival_time = np.random.exponential(l_arr[pt])
		new_arr = Patient(new_arrival_time + t[sim], t[sim], pt, 'arrival', 
			np.random.exponential(l_aban[pt]), np.random.exponential(w_mu[pt]))
		# Add new arrival to all queues
		for i in range(pm):

			heapq.heappush(A_Events[i], new_arr)
			a_queue_count[i][pt] += 1

# --------------------------------- Departure Event for non-continuous policies ---------------------------------
def departure_event(event, sim):

	global Events
	global A_Events
	global queue
	global w_mu
	global w_std
	global ward_capac
	global t
	global n_free
	global Treated
	global queue_length
	global ward_assignment
	global ward_nurse_def
	global k 
	global h_cost

	pt = event.get_pt()

	# Check if patients are in queue, checks if ward capacity allows for another patient
	if queue[sim][pt] and ward_capac[sim][pt] >= 0 and n_free[sim] >= 0:
		

		max_pt = pt

		# Get next patient and remove from Event list (abandoner) to change to a ward patient
		next_patient = queue[sim][max_pt].pop(0)
		# Events[sim].remove(next_patient)
		queue_length[sim][max_pt] -= 1
		
		# Change patient type to ward patient from abandoner
		next_patient.set_time(t[sim] + next_patient.get_serv())
		next_patient.set_location('ward')
		
		# Heapify the event list and replace the last departure with the new patient
		heapq.heappop(Events[sim])
		heapq.heappush(Events[sim], next_patient)

		ward_capac[sim][pt] += 1
		ward_capac[sim][max_pt] -= 1

	# Case when nurse finishes treating a patient and moves to her newly assigned ward
	elif ward_capac[sim][pt] < 0 and n_free[sim] >= 0 and ward_assignment[sim]:
		heapq.heappop(Events[sim])
		n_free[sim] += 1
		ward_capac[sim][pt] += 1

		w = random.choice(ward_assignment[sim])
		ward_nurse_def[sim][w] += 1
		if ward_nurse_def[sim][w] == 0:
			ward_assignment[sim].remove(w)

		if queue_length[sim][w] > 0:
			next_patient = queue[sim][w].pop(0)
			queue_length[sim][w] -= 1

			# Change patient type to ward patient from abandoner
			next_patient.set_time(t[sim] + next_patient.get_serv())
			next_patient.set_location('ward')
			
			# Heapify the event list and replace the last departure with the new patient
			heapq.heappush(Events[sim], next_patient)

			n_free[sim] -= 1
			ward_capac[sim][w] -= 1

	else:
		heapq.heappop(Events[sim])
		n_free[sim] += 1
		ward_capac[sim][pt] += 1

	Treated[sim] += 1

# --------------------------------- Abandonment Event for non-continuous policies ---------------------------------
def aban_event(event, sim):
	global Abandonment_Count
	global Events
	global A_Events
	global queue
	global queue_length

	pt = event.get_pt()


	abandoner = heapq.heappop(Events[sim])

	if abandoner in queue[sim][pt]:
		queue[sim][pt].remove(abandoner)
		Abandonment_Count[sim] += 1
		queue_length[sim][pt] -= 1



def simulation(T, N, lbda, mu, std, theta, tau, classes, hcost, q_cap, s_alloc, par_sim, rb, cont, preemption, s_t, p_t, a_t):
	global t, k, pm, r_time_arr, queue, l_arr, l_aban, w_mu, w_std, n_free, ward_alloc, ward_capac, last_arrival
	global Events, A_Events, capac, q_capac, a_queue_count, queue_length, ward_assignment, ward_nurse_def, h_cost, preempt
	

	# ----------------- Environment Variables -----------------
	l_arr = lbda 		# Arrival Rates
	l_aban = theta 		# Abandonment Rates
	w_mu = mu 			# Service Rate
	w_std = std 		# STD of Service Rate
	q_capac = q_cap 	# Capacity of each queue 
	r_time = tau 		# Shift Length
	Time = T 			# Total Simulation Run-Time
	pm = par_sim		# Total parallel simulations
	k = classes
	h_cost = hcost 
	preempt = preemption
	# ----------------- Simulation Variables -----------------
	# Variables Keeping track of states, queues, etc.
	
	# Initiate arrivals for first set of patients, same for all simulations
	e = [Patient(np.random.exponential(lbda[i]), t, i, 'arrival', 
		np.random.exponential(l_aban[i]), np.random.exponential(w_mu[i])) for i in range(0,k)]
	
	# Creating existing patients for the system
	existing = []
	for ind, e_p in enumerate(p_t):
		existing.append(Patient(0, 0, e_p, 'arrival', a_t[ind], s_t[ind]))

	# Sort patients in order of arrivals
	heapq.heapify(e)
	
	statistics = []

	for i in range(par_sim):
		t.append(0)				# Initiate Current time for each simulation, initial is 0
		queue.append([])  		# Initiate Queues for each simulation, initial is empty
		Events.append([]) 		# Initiate Event List for each simulation
		A_Events.append([])		# Initiate Arrival List
		n_free.append(N)		# Initiate Nurses free for each simulation, initial is N

		ward_alloc.append([]) 		# Initiate arrays for the allocation of nurses to each ward
		ward_capac.append([])		# Current Capacity of Each ward
		ward_nurse_def.append([])	# Keeps deficit of nurses during rebalance
		ward_assignment.append([])	# Keeps track of wards that still need nurses assigned during rebalance

		last_arrival.append([])		# Initiate arrays for each sim

		Balk_Count.append(0)			# Initiate Balk count for each sim
		Arrival_Count.append(0)			# Initiate Arrival count for each sim
		Abandonment_Count.append(0) 	# Initiate Abandonment count for each sim
		Treated.append(0)				# Initiate Treated patient count for each sim
		holding_cost.append(0)			# Initiate holding cost
		time_server_occupied.append(0)  # Initiate time occupied
		weighted_ward.append([]) 
		weighted_queue.append([]) 


		a_queue_count.append([]) 	# Initiate array for arrival queue count
		queue_length.append([]) 	# Initiate array for queue length

		r_time_arr.append(r_time) 	# Initiate first rebalance time

		statistics.append([]) 		# Initiate array recording arrivals

		headerX = []
		headerQ = []

		#Populating arrays
		for j in range(k):

			a_queue_count[i].append(1)			# Counts number of arrivals for a type of patient in the simulation
			queue_length[i].append(0)

			ward_alloc[i].append(s_alloc[i][j])	# Assigning the initial allocation for each ward
			ward_capac[i].append(s_alloc[i][j])	# Assigning the initial capacity free for each ward
			ward_nurse_def[i].append(0)			# Initiates nurse deficits in the ward to 0 at start

			last_arrival[i].append(0)			# Place holder for last arrival

			A_Events[i].append(e[j])  			# Assigning the initial patients to the events of each simulation

			queue[i].append([])

			weighted_queue[i].append(0)
			weighted_ward[i].append(0)			

			headerX.append('Ward_Count_' + str(j))
			headerQ.append('Queue_Count_' + str(j))

		for e_pat in existing:
			Events[i].append(e_pat)

		header = ['Time']
		statistics[i].append(header + headerX + headerQ)

	# ----------------- Simulation Start -----------------
	while(min(t) < T):
		for curr_sim in range(pm):
			
			#Check to see if we should continue running the current simulation
			if t[curr_sim] > T:
				continue
			row = []
			hc = 0
			servers_occupied = N - n_free[curr_sim]
			t_prev = 0

			# Calculating holding cost for next time period
			for wt in range(k):
				hc += queue_length[curr_sim][wt]*hcost[wt]


			# Check what the next event is
			if Events[curr_sim]:
				curr_event_E = Events[curr_sim][0]
				curr_event_A = A_Events[curr_sim][0]
				if curr_event_E.get_time() < curr_event_A.get_time():
					curr_event = curr_event_E
				else:
					curr_event = curr_event_A
			else:
				curr_event = A_Events[curr_sim][0]


			# Main body of simulation
			if cont[curr_sim] == 1:
				t_prev = t[curr_sim]
				t[curr_sim] = curr_event.get_time()
				
				if curr_event.get_location() == 'arrival':
					# print 'arrival:: Nurses: ' + str(n_free) + ' Queue Length:' + str(queue_length) + ' ' + 'Ward capac:' + str(ward_capac) + 'Ward alloc:' + str(ward_alloc) + ' ' + str(curr_sim) + ' ' + str(t[curr_sim])
					arrival_event_cont(curr_event, curr_sim)
					
				elif curr_event.get_location() == 'ward':
					# print 'departure:: Nurses: ' + str(n_free) + ' Queue Length:' + str(queue_length) + ' ' + 'Ward capac:' + str(ward_capac) + 'Ward alloc:' + str(ward_alloc) + ' ' + str(curr_sim) + ' ' + str(t[curr_sim])
					departure_event_cont(curr_event, curr_sim)
					
				elif curr_event.get_location() == 'abandonment':
					# print 'arrival:: Nurses: ' + str(n_free) + ' Queue Length:' + str(queue_length) + ' ' + 'Ward capac:' + str(ward_capac) + 'Ward alloc:' + str(ward_alloc) + ' ' + str(curr_sim) + ' ' + str(t[curr_sim])
					aban_event(curr_event, curr_sim)
			else:
				
				# Rebalance case
				if r_time_arr[curr_sim] < curr_event.get_time() and rb[curr_sim] == 1:
					t_prev = t[curr_sim]
					t[curr_sim] = r_time_arr[curr_sim]
					rebalance(N, curr_sim)
					# print 'rebalance ' + str(ward_nurse_def[curr_sim])
					r_time_arr[curr_sim] = r_time_arr[curr_sim] + r_time


				# No Rebalance case
				elif r_time_arr[curr_sim] < curr_event.get_time() and rb[curr_sim] == 0:
					t_prev = t[curr_sim]
					t[curr_sim] = r_time_arr[curr_sim]
					r_time_arr[curr_sim] = r_time_arr[curr_sim] + r_time

				# Continue
				else:
					#print 'else'
					t_prev = t[curr_sim]
					t[curr_sim] = curr_event.get_time()
					
					if curr_event.get_location() == 'arrival':
						# print 'arrival:: Nurses: ' + str(n_free) + ' Queue Length:' + str(queue_length) + ' ' + 'Ward capac:' + str(ward_capac) + ' Ward alloc:' + str(ward_alloc) + ' ' + str(curr_sim) + ' ' + str(t[curr_sim])
						arrival_event(curr_event, curr_sim)
						
					elif curr_event.get_location() == 'ward':
						# print 'arrival:: Nurses: ' + str(n_free) + ' Queue Length:' + str(queue_length) + ' ' + 'Ward capac:' + str(ward_capac) + ' Ward alloc:' + str(ward_alloc) + ' ' + str(curr_sim) + ' ' + str(t[curr_sim])
						departure_event(curr_event, curr_sim)
						
					elif curr_event.get_location() == 'abandonment':
						# print 'arrival:: Nurses: ' + str(n_free) + ' Queue Length:' + str(queue_length) + ' ' + 'Ward capac:' + str(ward_capac) + ' Ward alloc:' + str(ward_alloc) + ' ' + str(curr_sim) + ' ' + str(t[curr_sim])
						aban_event(curr_event, curr_sim)
			
			holding_cost[curr_sim] += (t[curr_sim] - t_prev)*hc
			time_server_occupied[curr_sim] +=  (t[curr_sim] - t_prev)*servers_occupied

			row.append(t[curr_sim])

			w_val = 0
			for x in range(k):
				w_val = ward_alloc[curr_sim][x]-ward_capac[curr_sim][x]+queue_length[curr_sim][x]
				row.append(w_val)
				weighted_ward[curr_sim][x] += (t[curr_sim] - t_prev)*w_val
			for x in range(k):
				w_val = queue_length[curr_sim][x]
				row.append(w_val)
				weighted_queue[curr_sim][x] += (t[curr_sim] - t_prev)*w_val
			statistics[curr_sim].append(row)
	# print len(statistics[0])
	return statistics

# Function to write table results to a file, used in the main function to write the return values to a csv file
def writeLog(fil, table):
	c1 = csv.writer(fil)

	for val in table:
		# print val
		c1.writerow(val)


Trials = 10
Results = []



# ================ Input Variables ================
Total_Time = 400
lbda_out = [1.0/15.0, 1.0/15.0]
mu_out = [1.0/8.0, 1.0/8.0]
std_out = [1, 1]
theta_out = [10000, 10000]
tau_out = .5
k_out = 2
Nurses = 4
hcost_out = [1,2]
q_cap_out = [float('inf'), float('inf')]

# Parallel simulation variables
s_alloc_out = [[2,2], [4,4]]
rebalance1 = [1, 1]
cont_out = [0, 1]
preemption_out = [0, 1]

stats = []

# ---------------- Pre-existing Patients in ward ----------------
service_times = []
patient_type =  []
abandonment_times = []
patient_type_count = [5,5]

for ind, cnt in enumerate(patient_type_count):
	service_times = service_times + [np.random.exponential(mu_out[ind]) for x in range(cnt)]
	patient_type = patient_type + [ind for x in range(cnt)]
	abandonment_times += [10000 for x in range(cnt)]


	# print service_times
	# print patient_type
	# print abandonment_times

	# service_times = [.15, .05, .21, .04]
	# patient_type =  [1, 0, 1, 0]
	# abandonment_times = [10000, 10000, 10000, 10000]

	
for tests in range(Trials):
	

	stats = simulation(Total_Time, Nurses, lbda_out, mu_out, std_out, theta_out, tau_out, k_out, hcost_out, q_cap_out, s_alloc_out, 2, rebalance1, cont_out, preemption_out, service_times, patient_type, abandonment_times)


	resetvar()

	# print Treated
	# print Arrival_Count
	# print Abandonment_Count
	# print Balk_Count
	# print holding_cost
	# print time_server_occupied
	# print Total_Time*Nurses

	print "Done with Trial " + str(tests)

# ================ File Writing ================

fil0 = open(os.getcwd() + "/Sim_Rebalance_20.csv","wb")
fil1 = open(os.getcwd() + "/Sim_No_Rebalance_20.csv","wb")

writeLog(fil0, stats[0])
writeLog(fil1, stats[1])

# print Counter_Record

print '\n'

for k in range(0,2):
	dataset_arr = []
	dataset_hc = []
	dataset_st = []
	dataset_wq = []
	dataset_ww = []
	for i in range(k_out):
		dataset_wq.append([])
		dataset_ww.append([])

	print "Simulation " + str(k) 
	for i in range(Trials):
		dataset_arr.append(Counter_Record[i][0][k])	
		dataset_hc.append(Counter_Record[i][1][k])
		dataset_st.append(Counter_Record[i][2][k])
		for j in range(k_out):
			dataset_wq[j].append(Counter_Record[i][3][k][j]/Counter_Record[i][5][k])
			dataset_ww[j].append(Counter_Record[i][4][k][j]/Counter_Record[i][5][k])

	print "Arrivals CI: " + str(mean_confidence_interval(dataset_arr))
	print "Holding Cost CI: " + str(mean_confidence_interval(dataset_hc))
	print "Server Time CI: " + str(mean_confidence_interval(dataset_st))
	for i in range(k_out):
		print "Queue length for ward CI " + str(i) + ": " + str(mean_confidence_interval(dataset_wq[i]))
		print "Headcount for ward CI" + str(i) + ": " + str(mean_confidence_interval(dataset_ww[i]))
	print '\n'




