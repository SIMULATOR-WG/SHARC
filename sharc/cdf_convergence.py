import pickle
import os
import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

#output_dir_vec = ["output_2018-03-13_01", "output_2018-03-13_02"]
output_dir_vec = ["output_2018-03-14_01", "output_2018-03-14_02", "output_2018-03-14_03"]

k = len(output_dir_vec)

# significance_level = 0.05, i.e. probability of false negative (H0 = distributions are the same)
t_vec = [1.96, 1.945, 1.915, 1.894]
if 2 <= k <= 5:
    t = t_vec[k-2]
else:
    b0 = 1.645
    b1 = .678
    b2 = -.362
    t = b0 + b1/np.sqrt(k-1) + b2/(k-1)

range_max = 20000

# n = number of snapshots
# k = number of samples
n_vec = np.concatenate( (np.arange(50, 500, 50), np.arange(500, range_max + 1, 500)))
AD_metric_vec = np.array([])
KS_metric_vec = np.array([])
DB_metric_vec = np.array([])

for n in n_vec:

    msg = "number of drops = " + str(n)
    print(msg)

    all_samples = np.empty([len(output_dir_vec), n])
    index = 0

    for output_dir in output_dir_vec:
        filename = 'ordered_sample_' + str(n) + ".dat"
        filename = os.path.join(output_dir, filename)

        with open(filename, 'rb') as fp:
            sample = np.array(pickle.load(fp))
            all_samples[index,:] = sample

        index += 1

    ####################################################################
    # Anderson-Darling Test (considering all samples have the same size)

    # # build pooled ordered samples
    # ordered_samples = np.sort(all_samples.flatten())
    #
    # # calculate A_kN^2
    # N = ordered_samples.size
    # A_kn2 = 0
    #
    # for i in range(k):
    #     for j in range(1, N):
    #         Mij = np.searchsorted(all_samples[i,:], ordered_samples[j-1], side="right")
    #         A_kn2 += (N*Mij-j*n)**2/ (j *(N-j)) / n / N
    #
    # H = k / n
    # h = np.sum(1/np.arange(1, N))
    #
    # g = 0
    # for i in range(1, N-1):
    #     aux = np.arange(i+1, N)
    #     g+= np.sum(1/aux/(N-i))
    #
    # a = (4*g-6)*(k-1) + (10-6*g)*H
    # b = (2*g-4)*k**2 + 8*h*k + (2*g-14*h-4)*H - 8*h + 4*g - 6
    # c = (6*h+2*g-2)*k**2 + (4*h-4*g+6)*k + (2*h-6)*H + 4*h
    # d = (2*h+6)*k**2 - 4*h*k
    #
    # sigma2 = (a*N**3 + b*N**2 + c*N + d)/((N-1)*(N-2)*(N-3))
    #
    # Tkn = (A_kn2 - (k-1))/np.sqrt(sigma2)
    #
    # AD_metric_vec = np.append(AD_metric_vec, Tkn)


    ####################################################################
    # Smirnov-Kolmogorov Test (considering all samples have the same size)
    KS_metric = 0
    number_of_pairs = 0

    cdf = (np.arange(n)+1)/n

    min_x = np.min(all_samples, 1)
    min_x = np.max(min_x)

    max_x = np.max(all_samples, 1)
    max_x = np.min(max_x)

    x = np.linspace(min_x, max_x, n)

    for first_sample_index in range(k):
        for second_sample_index in range (first_sample_index + 1, k):
            number_of_pairs += 1
            f1 = interp1d(all_samples[first_sample_index,:], cdf, kind='slinear')
            f2 = interp1d(all_samples[second_sample_index,:], cdf, kind='slinear')

            cdf1 = f1(x)
            cdf2 = f2(x)

            KS_metric += np.max(abs(cdf1-cdf2))

    KS_metric = KS_metric / number_of_pairs
    KS_metric_vec = np.append(KS_metric_vec, KS_metric)

    ####################################################################
    # Bhattacharyya distance (considering all samples have the same size)
    DB_metric = 0
    number_of_pairs = 0

    hist_all, bin_edges = np.histogram(all_samples)
    pmf_all = hist_all / (n*k)

    for index in range(k):
        hist, dummy = np.histogram(all_samples[index,:], bin_edges)
        pmf = hist / n
        BC = np.sum(np.sqrt(pmf_all * pmf))
        DB_metric += -np.log(BC)

    DB_metric = DB_metric / k
    DB_metric_vec = np.append(DB_metric_vec, DB_metric)

#
# plt.plot(n_vec, AD_metric_vec, n_vec, np.ones(n_vec.size) * t)
# plt.xlabel("number of drops")
# plt.ylabel("Anderson-Darling statistic")

plt.plot(n_vec, KS_metric_vec)
plt.xlabel("number of drops")
plt.ylabel("Kolmogorov-Smirnov statistic")


plt.figure()
plt.plot(n_vec, DB_metric_vec)
plt.xlabel("number of drops")
plt.ylabel("average Bhattacharyya distance")

plt.show()
