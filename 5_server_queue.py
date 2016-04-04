import numpy as np 
import heapq
import os
import time
import random 
import csv

# Global Variables for easier use in the simulation.

pm = 0				# Number of paral
k = 0
t = []
r_time_arr = []
queue = []
l_arr = []
l_aban = []
w_mu = []
w_std = []
n_free = []

q_capac = []
Events = []

Arrival_Count = []
Abandonment_Count = []
Balk_Count = []
Treated = []
holding_cost = []
time_server_occupied = []

a_queue_count = []

queue_length = []
ward_alloc = []
ward_capac = []
ward_nurse_def = []
ward_assignment = []

arrival_queue = []

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



def rebalance(N, sim):

	global queue_length

	global ward_capac
	global ward_alloc
	global ward_assignment
	global ward_nurse_def
	global Events
	global n_free 
	global pm
	global k

	print queue_length[sim]
	print ward_capac[sim]
	print n_free

	old_alloc = [0]*k
	total = max(sum(queue_length[sim]), k)
	s = 0
	for i in range(k-1):
		old_alloc[i] = ward_alloc[sim][i]
		ward_alloc[sim][i] = queue_length[sim][i]*(N-k)/total + 1
		ward_capac[sim][i] = ward_capac[sim][i] + ward_alloc[sim][i] - old_alloc[i]
		s += queue_length[sim][i]*(N-k)/total


	old_alloc[k-1] = ward_alloc[sim][k-1]
	ward_alloc[sim][k-1] = (N - k) - s + 1
	ward_capac[sim][k-1] = ward_capac[sim][k-1] + ward_alloc[sim][k-1] - old_alloc[k-1]
	
	print 'phase i'
	print ward_capac[sim]
	print ward_alloc[sim]
	print n_free[sim]
	
	for i in range(k):
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
	for i in range(k):
		if n_free[sim] == 0 and ward_capac[sim][i] > 0:
			ward_nurse_def[sim][i] = -ward_capac[sim][i]
			ward_assignment[sim].append(i)
	
	print 'phase ii'
	print ward_capac[sim]
	print ward_alloc[sim]
	print n_free[sim]

def arrival_event(event, sim):
	global Arrival_Count
	global Balk_Count
	global ward_capac
	global w_mu
	global w_std
	global t
	global l_aban
	global Events
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
		# heapq.heappush(Events[sim], new_queue_arr)
		queue[sim][pt].append(new_queue_arr)
		queue_length[sim][pt] += 1

	#Case where queue is full and patient leaves system
	else:
		Balk_Count[sim] += 1

	# Remove Patient Arrival Event and reduce count
	heapq.heappop(Events[sim])
	a_queue_count[sim][pt] -= 1

	# print a_queue_count[sim]
	#Create new arrival and replace old arrival
	if a_queue_count[sim][pt] == 0:

		# Generate new arrival time for patient and create new patient
		new_arrival_time = np.random.exponential(l_arr[pt])
		new_arr = Patient(new_arrival_time + t[sim], t[sim], pt, 'arrival', 
			np.random.exponential(l_aban[pt]), np.random.exponential(w_mu[pt]))
		# Add new arrival to all queues
		for i in range(pm):

			heapq.heappush(Events[i], new_arr)
			a_queue_count[i][pt] += 1


def departure_event(event, sim):

	global Events
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

	pt = event.get_pt()

	# Check if patients are in queue, checks if ward capacity allows for another patient
	if queue[sim][pt] and ward_capac[sim][pt] >= 0 and n_free[sim] >= 0:
		# Get next patient and remove from Event list (abandoner) to change to a ward patient
		next_patient = queue[sim][pt].pop(0)
		# Events[sim].remove(next_patient)
		queue_length[sim][pt] -= 1
		
		# Change patient type to ward patient from abandoner
		next_patient.set_time(t[sim] + next_patient.get_serv())
		next_patient.set_location('ward')
		
		# Heapify the event list and replace the last departure with the new patient
		heapq.heappop(Events[sim])
		heapq.heappush(Events[sim], next_patient)
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

def aban_event(event, sim):
	global Abandonment_Count
	global Events
	global queue
	global queue_length

	pt = event.get_pt()


	abandoner = heapq.heappop(Events[sim])

	if abandoner in queue[sim][pt]:
		queue[sim][pt].remove(abandoner)
		Abandonment_Count[sim] += 1
		queue_length[sim][pt] -= 1



def simulation(T, N, lbda, mu, std, theta, tau, classes, hcost, q_cap, s_alloc, par_sim, rb):
	global t, k, pm, r_time_arr, queue, l_arr, l_aban, w_mu, w_std, n_free, ward_alloc, ward_capac, Events, capac, q_capac, a_queue_count, queue_length, ward_assignment, ward_nurse_def
	

	# ----------------- Environment Variables -----------------
	l_arr = lbda 		# Arrival Rates
	l_aban = theta 		# Abandonment Rates
	w_mu = mu 			# Service Rate
	w_std = std 		# STD of Service Rate
	q_capac = q_cap 	# Capacity of each queue 
	r_time = tau 		# Shift Length
	Time = T 			# Total Simulation Run-Time
	pm = par_sim			# Total parallel simulations
	k = classes 
	# ----------------- Simulation Variables -----------------
	# Variables Keeping track of states, queues, etc.
	
	# Initiate arrivals for first set of patients, same for all simulations

	e = [Patient(np.random.exponential(lbda[i]), t, i, 'arrival', 
		np.random.exponential(l_aban[i]), np.random.exponential(w_mu[i])) for i in range(0,k)]

	# Sort patients in order of arrivals
	heapq.heapify(e)
	
	statistics = []

	for i in range(par_sim):
		t.append(0)				# Current time for each simulation, initial is 0
		queue.append([])  		# Queues for each simulation, initial is empty
		Events.append([]) 		# Event List for each simulation
		n_free.append(N)		# Nurses free for each simulation, initial is N

		ward_alloc.append([]) 		# The allocation of nurses to each ward
		ward_capac.append([])		# Current Capacity of Each ward
		ward_nurse_def.append([])	# Keeps deficit of nurses during rebalance
		ward_assignment.append([])	# Keeps track of wards that still need nurses assigned during rebalance

		Balk_Count.append(0)			# Balk count for each sim
		Arrival_Count.append(0)			# Arrival count for each sim
		Abandonment_Count.append(0) 	# Abandonment count for each sim
		Treated.append(0)				# Treated patient count for each sim
		holding_cost.append(0)
		time_server_occupied.append(0)

		a_queue_count.append([])
		queue_length.append([])

		r_time_arr.append(r_time)

		statistics.append([])

		headerX = []
		headerQ = []
		for j in range(k):

			a_queue_count[i].append(1)			# Counts number of arrivals for a type of patient in the simulation
			queue_length[i].append(0)

			ward_alloc[i].append(s_alloc[i][j])	# Assigning the initial allocation for each ward
			ward_capac[i].append(s_alloc[i][j])	# Assigning the initial capacity free for each ward
			ward_nurse_def[i].append(0)			# Initiates nurse deficits in the ward to 0 at start

			Events[i].append(e[j])  			# Assigning the initial patients to the events of each simulation

			queue[i].append([])

			headerX.append('Ward_Count_' + str(j))
			headerQ.append('Queue_Count_' + str(j))

		header = ['Time']
		statistics[i].append(header + headerX + headerQ)

	# ----------------- Simulation Start -----------------
	while(min(t) < T):
		for curr_sim in range(pm):
			if t[curr_sim] > T:
				continue
			row = []
			hc = 0
			servers_occupied = N - n_free[curr_sim]
			t_prev = 0
			for wt in range(k):
				hc += queue_length[curr_sim][wt]*hcost[wt]
			

			curr_event = Events[curr_sim][0]
			if r_time_arr[curr_sim] < curr_event.get_time() and rb[curr_sim] == 1:
				t_prev = t[curr_sim]
				t[curr_sim] = r_time_arr[curr_sim]
				rebalance(N, curr_sim)
				print 'rebalance ' + str(ward_nurse_def[curr_sim])
				r_time_arr[curr_sim] = r_time_arr[curr_sim] + r_time


			elif r_time_arr[curr_sim] < curr_event.get_time() and rb[curr_sim] == 0:
				t_prev = t[curr_sim]
				t[curr_sim] = r_time_arr[curr_sim]
				r_time_arr[curr_sim] = r_time_arr[curr_sim] + r_time


			else:
				#print 'else'
				t_prev = t[curr_sim]
				t[curr_sim] = curr_event.get_time()
				
				if curr_event.get_location() == 'arrival':
					print 'arrival ' + str(n_free) + ' ' + str(queue_length) + ' ' + str(curr_sim) 
					arrival_event(curr_event, curr_sim)
					
				elif curr_event.get_location() == 'ward':
					print 'ward ' + str(n_free) + ' ' + str(queue_length) + ' ' + str(curr_sim) 
					departure_event(curr_event, curr_sim)
					
				elif curr_event.get_location() == 'abandonment':
					print 'abandonment ' + str(n_free) + ' ' + str(queue_length) + ' ' + str(curr_sim) 
					aban_event(curr_event, curr_sim)
			
			holding_cost[curr_sim] += (t[curr_sim] - t_prev)*hc
			time_server_occupied[curr_sim] +=  (t[curr_sim] - t_prev)*servers_occupied

			row.append(t[curr_sim])
			for x in range(k):
				row.append(ward_alloc[curr_sim][x]-ward_capac[curr_sim][x]+queue_length[curr_sim][x])	
			for x in range(k):
				row.append(queue_length[curr_sim][x])
			statistics[curr_sim].append(row)
	print len(statistics[0])
	return statistics
# Function to write table results to a file, used in the main function to write the return values to a csv file
def writeLog(fil, table):
	c1 = csv.writer(fil)

	for val in table:
		# print val
		c1.writerow(val)

Total_Time = 400
lbda_out = [1.0/15.0, 1.0/15.0]
mu_out = [1.0/8.0, 1.0/8.0]
std_out = [1, 1]
theta_out = [10000, 10000]
tau_out = .5
k_out = 2
hcost_out = [1,1]
q_cap_out = [float('inf'), float('inf')]
s_alloc_out = [[2,2], [2,2]]
rebalance1 = [1, 0]
stats = []
Nurses = 4
stats = simulation(Total_Time, Nurses, lbda_out, mu_out, std_out, theta_out, tau_out, k_out, hcost_out, q_cap_out, s_alloc_out, 2, rebalance1)

fil0 = open(os.getcwd() + "/Sim_Rebalance_20.csv","wb")
fil1 = open(os.getcwd() + "/Sim_No_Rebalance_20.csv","wb")

writeLog(fil0, stats[0])
writeLog(fil1, stats[1])


print Treated
print Arrival_Count
print Abandonment_Count
print Balk_Count
print holding_cost
print time_server_occupied
print Total_Time*Nurses
