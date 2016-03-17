import numpy as np 
import heapq
import os
import time

# Global Variables for easier use in the simulation.

t = 0
r_time = 12
queue = []
l_arr = []
l_aban = []
w_mu = []
w_std = []
n_free = 0
ward_alloc = []
ward_capac = []
q_capac = []
Events = []

Arrival_Count = 0
Abandonment_Count = 0
Balk_Count = 0
Treated = 0

class Patient:

	def __init__(self, t, t_arr, pt, location):
		self.t = t
		self.location = location
		self.ward = 0
		self.t_arr = t_arr
		self.pt = pt

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



def rebalance(q1,q2,q3,q4,q5, N):

	n1 = N/5
	n2 = N/5
	n3 = N/5
	n4 = N/5
	n5 = N/5

	n_array = [n1,n2,n3,n4,n5]

	return n_array

def arrival_event(event):
	global Arrival_Count
	global BalkCount
	global ward_capac
	global w_mu
	global w_std
	global t
	global l_aban
	global Events
	global l_arr
	global queue 
	global n_free

	Arrival_Count += 1

	pt = event.get_pt()
	
	if ward_capac[pt] > 0 and n_free > 0:
		new_ward_arr = Patient(np.random.lognormal(w_mu[pt], w_std[pt]) + t, t, pt, 'ward')
		heapq.heappush(Events, new_ward_arr)
		ward_capac[pt] = ward_capac[pt] - 1
		n_free = n_free - 1
	elif len(queue[pt]) < q_capac[pt]:
		new_queue_arr = Patient(np.random.exponential(l_aban[pt]) + t, t, pt, 'abandonment')
		heapq.heappush(Events, new_queue_arr)
		queue[pt].append(new_queue_arr)
	else:
		Balk_Count += 1

	new_arr = Patient(np.random.exponential(l_arr[pt]) + t, t, pt, 'arrival')
	heapq.heappushpop(Events, new_arr)
	
	

def departure_event(event):

	global Events
	global queue
	global w_mu
	global w_std
	global ward_capac
	global t
	global n_free
	global Treated

	pt = event.get_pt()

	if queue[pt] and ward_capac > 0 and n_free > 0:
		next_patient = queue[pt].pop(0)
		Events.remove(next_patient)
		heapq.heapify(Events)
		heapq.heappushpop(Events, next_patient)
	else:
		heapq.heappop(Events)
		n_free += 1
		ward_capac[pt] += 1

	Treated += 1

def aban_event(event):
	global Abandonment_Count
	global Events
	global queue

	pt = event.get_pt()


	abandoner = heapq.heappop(Events)

	if abandoner in queue[pt]:
		queue[pt].remove(abandoner)
		Abandonment_Count += 1

def simulation(T, N, l0, l1, l2, l3, l4, mu0, mu1, mu2, mu3, mu4, std0, std1, std2, std3, std4, a0, a1, a2, a3, a4):
	global t, r_time, queue, l_arr, l_aban, w_mu, w_std, n_free, ward_alloc, ward_capac, Events, q_capac

	t = 0
	r_time = 12
	Time = T
	
	# queues
	queue = [[], [], [], [], []]
	q_capac = [6, 6, 6, 6, 6]

	l_arr = [l0, l1, l2, l3, l4]
	l_aban = [a0, a1, a2, a3, a4]
	w_mu = [mu0, mu1, mu2, mu3, mu4]
	w_std = [std0, std1, std2, std3, std4]

	n_free = N

	# nurses in each ward / array or dict?
	ward_alloc = [4, 4, 4, 4, 4]
	ward_capac = [4, 4, 4, 4, 4]

	Events = []

	Events.append(Patient(np.random.exponential(l0), t, 0, 'arrival'))
	Events.append(Patient(np.random.exponential(l1), t, 1, 'arrival'))
	Events.append(Patient(np.random.exponential(l2), t, 2, 'arrival'))
	Events.append(Patient(np.random.exponential(l3), t, 3, 'arrival'))
	Events.append(Patient(np.random.exponential(l4), t, 4, 'arrival'))

	heapq.heapify(Events)

	while(t < T):
		curr_event = Events[0]

		if r_time < curr_event.get_time():
			t = r_time
			# rebalance()
			r_time = r_time + 12
		else:
			#print 'else'
			t = curr_event.get_time()
			
			if curr_event.get_location() == 'arrival':
				arrival_event(curr_event)
				print 'arrival'
			elif curr_event.get_location() == 'ward':
				departure_event(curr_event)
				print 'ward'
			elif curr_event.get_location() == 'abandonment':
				aban_event(curr_event)
				print 'abandonment'



simulation(100, 20, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
print Treated
print Arrival_Count
print Abandonment_Count
print Balk_Count
