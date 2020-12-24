import numpy as np
import matplotlib.pyplot as plt
# See https://www.researchgate.net/publication/226675412_A_Probabilistic_Interpretation_of_Precision_Recall_and_F-Score_with_Implication_for_Evaluation

def posterior_distribution_f1_score(conf_mat):
    tp, fp, fn  = conf_mat[1][1], conf_mat[0][1], conf_mat[1][0]
    lbda = 1
    h =1
    u = np.random.gamma(size=1000000, scale=2*h, shape=tp + lbda)
    v = np.random.gamma(size=1000000, scale=h, shape=fp + fn + 2*lbda)
    out = u/(u + v)
    print("Mean: {}, Stddev: {}".format(np.mean(out), np.std(out)))
    _ = plt.hist(out, bins='auto')