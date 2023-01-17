import numpy as np
import matplotlib.pyplot as plt


#QE = 1#.8
Nr = 1.4
P = 1.
#t = .1
D = 0#.5

def SNR(t,QE,Nr,D):
	return QE*P*t/(np.sqrt(QE*P*t+D*t+Nr**2.))

t_lst = np.linspace(0.1,100,1000)
SNR_lst_Nr0D0 = []
SNR_lst_Nr14D1 = []
SNR_lst_Nr7D5 = []
#SNR_lst_ = []

for i in range(len(t_lst)):
    SNR_lst_Nr0D0.append(SNR(t_lst[i],1,0,0))
    SNR_lst_Nr14D1.append(SNR(t_lst[i],.8,1.4,.1/2.))
    SNR_lst_Nr7D5.append(SNR(t_lst[i],.8,.7,.5/2.))

fig, ax = plt.subplots()

ax.set(title='SNR vs photons/pixel/second')


# Plotting the graph with Log ticks at x and y axis using loglog
ax.loglog(t_lst, SNR_lst_Nr0D0, '--r', linewidth=2, label=r'SNR QE=1 N$_r$ =0.0 D=0.0')
ax.loglog(t_lst, SNR_lst_Nr14D1, '--g', linewidth=2, label=r'SNR QE=0.8 N$_r$ =1.4 D=0.1')
ax.loglog(t_lst, SNR_lst_Nr7D5, '--b', linewidth=2, label=r'SNR QE=0.8 N$_r$ =0.7 D=0.5')

ax.set_title('SNR vs. photon count (rate: 1/s)', fontsize=15)
ax.set_xlabel('Exposure time [s]', fontsize=13)
ax.set_ylabel('SNR', fontsize=13)
ax.legend()

plt.tight_layout()
plt.show()

