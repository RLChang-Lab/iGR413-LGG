################## GLMM validation of Threshold for IV DO data ##################
#   This code takes the fold change values from the IV data, plots a distribution, and identifies intersections that describe the seperation of the data in the distribution
#   Seeks to answer the question: what fold change splits the data set in 2 groups? 
#   Step 1: load in the fold change data in xl format
#   Step 2: Run plot_gmm_histogram
        #Fits 2 Gaussian distributions
        # find_gaussian_intersection identifies the intersection point and assigns a confidence score to the intersection
        #plots the distribution, the 2 distributions, and lists the threshold and confidenc escore 

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.mixture import GaussianMixture
from scipy.stats import norm

def find_gaussian_intersection(m1, s1, w1, m2, s2, w2):
    a = 1/(2*s1**2) - 1/(2*s2**2)
    b = m2/(s2**2) - m1/(s1**2)
    c = (m1**2)/(2*s1**2) - (m2**2)/(2*s2**2) + np.log((s2*w1)/(s1*w2))
    roots = np.roots([a, b, c])
    real_roots = roots[np.isreal(roots)].real
    real_roots = real_roots[(real_roots > min(m1, m2)) & (real_roots < max(m1, m2))]
    return real_roots

def plot_gmm_histogram(data, title, filename=None):
    data = data[~np.isnan(data)] 
    values = data.values.reshape(-1,1)
    
    if len(values) < 10:
        print(f"Skipping {title}, not enough data points.")
        return

    gmm = GaussianMixture(n_components=2, random_state=42)
    gmm.fit(values)

    means = gmm.means_.flatten()
    variances = gmm.covariances_.flatten()
    weights = gmm.weights_
    stds = np.sqrt(variances)
    
    order = np.argsort(means)
    means = means[order]
    variances = variances[order]
    weights = weights[order]
    stds = stds[order]

    intersections = find_gaussian_intersection(means[0], stds[0], weights[0],
                                               means[1], stds[1], weights[1])
    threshold = intersections[0] if len(intersections) > 0 else None

    x = np.linspace(min(values), max(values), 1000).flatten()
    pdf1 = weights[0] * norm.pdf(x, means[0], stds[0])
    pdf2 = weights[1] * norm.pdf(x, means[1], stds[1])
    overlap_area = np.trapz(np.minimum(pdf1, pdf2), x)
    confidence = 1 - overlap_area  
    plt.figure(figsize=(7,5))
    plt.hist(values.flatten(), bins=30, density=True, alpha=0.5, color='gray')
    plt.plot(x, pdf1, label='GMM Component 1')
    plt.plot(x, pdf2, label='GMM Component 2')
    plt.plot(x, pdf1 + pdf2, linestyle='--', label='GMM Combined')
    if threshold is not None:
        plt.axvline(threshold, color='red', linestyle=':', label=f'Threshold: {threshold:.3f}')
    plt.title(f"{title}\nSeparation confidence: {confidence:.2f} (1=perfect)")
    plt.xlabel('Fold Change')
    plt.ylabel('Density')
    plt.legend()
    if filename:
        plt.savefig(filename, bbox_inches='tight')
    plt.show()
    
    return threshold, confidence

df = pd.read_excel('iv_fc_paper.xlsx', index_col=0)
combined_data = df.melt(var_name='Treatment', value_name='Value')['Value']
plot_gmm_histogram(combined_data, "Combined DM Dropouts (DM52, DM24, DM16, DM13)", filename='combined_gmm_hist_final.svg')
