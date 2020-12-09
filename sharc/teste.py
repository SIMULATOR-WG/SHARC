import numpy as np

#d=20000
#f=2150

#loss = 20*np.log10(d) + 20*np.log10(f) - 27.55

#print(loss)
#a= (np.arange(.05, 5, .1) - .05)
#print(a)
emission_limits = -7 - 7 / 5 * (np.arange(.05, 5, .1) - .05)
            #emission_limits = -7 - 7 / 5 * (np.arange(.05, 5, .1) - .05)
print(emission_limits)
emission_limits = np.append(emission_limits, np.array([-14, -14]))
print(emission_limits)
