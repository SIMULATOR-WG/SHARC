import numpy as np # for numerical computing
import matplotlib.pyplot as plt # for plotting functions
from matplotlib import cm # colormap for color palette
from scipy.special import erfc
from sharc.support.modulation.modem import PSKModem,QAMModem,PAMModem,FSKModem
from sharc.support.modulation.channels import awgn
from sharc.support.modulation.errorRates import ser_awgn

#---------Input Fields------------------------
nSym = 10**6 # Number of symbols to transmit
EbN0dBs = np.arange(start=-4,stop = 24, step = 2) # Eb/N0 range in dB for simulation
mod_type = 'QAM' # Set 'PSK' or 'QAM' or 'PAM' or 'FSK'
#arrayOfM = [2,4,8,16,32] # array of M values to simulate
arrayOfM=[16, 64, 256] # uncomment this line if MOD_TYPE='QAM'
#arrayOfM=[2, 4] # uncomment this line if MOD_TYPE='QAM'
coherence = 'coherent' #'coherent'/'noncoherent'-only for FSK
modem_dict = {'psk': PSKModem,'qam':QAMModem,'pam':PAMModem,'fsk':FSKModem}
colors = plt.cm.jet(np.linspace(0,1,len(arrayOfM))) # colormap
fig, ax = plt.subplots(nrows=1,ncols = 1)
for i, M in enumerate(arrayOfM):
    # -----Initialization of various parameters----
    k = np.log2(M)
    EsN0dBs = 10 * np.log10(k) + EbN0dBs  # EsN0dB calculation
    SER_sim = np.zeros(len(EbN0dBs))  # simulated Symbol error rates
    inputSyms = np.random.randint(low=0, high=M, size=nSym)
    # uniform random symbols from 0 to M-1
    if mod_type.lower() == 'fsk':
        modem = modem_dict[mod_type.lower()](M)  # choose modem from dictionary
    else:  # for all other modulations
        modem = modem_dict[mod_type.lower()](M)  # choose modem from dictionary
    modulatedSyms = modem.modulate(inputSyms)  # modulate

    for j, EsN0dB in enumerate(EsN0dBs):
        receivedSyms = awgn(modulatedSyms, EsN0dB)  # add awgn noise
        if mod_type.lower() == 'fsk':  # demodulate (Refer Chapter 3)
            detectedSyms = modem.demodulate(receivedSyms,coherence)
        else:  # demodulate (Refer Chapter 3)
            detectedSyms = modem.demodulate(receivedSyms)
        SER_sim[j] = np.sum(detectedSyms != inputSyms) / nSym

    SER_theory = ser_awgn(EbN0dBs, mod_type, M, coherence)  # theory SER
    ax.semilogy(EbN0dBs, SER_sim, color= colors[i], marker='o', linestyle='', label='Sim ' + str(M) + '-' + mod_type.upper())
    ax.semilogy(EbN0dBs, SER_theory, color=colors[i], linestyle='-', label='Theory' + str(M) +'-'+mod_type.upper())

M = 4 #16 points in the constellation
pskModem = PSKModem(M) #create a 16-PSK modem object
pskModem.plotConstellation() #plot ideal constellation for this modem
ax.set_xlabel('Eb/N0(dB)')
ax.set_ylabel('SER ($P_s$)')
plt.ylim(10E-7, 0)
plt.xlim(-4, 24)
ax.set_title('Probability of Symbol Error for M-'+str(mod_type)+' over AWGN')
ax.legend()
plt.show()
