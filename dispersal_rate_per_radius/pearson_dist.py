import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import pearson3, norm

# Load data
df = pd.read_csv("empirical_dispersal_moments_grouped.tsv", sep="\t")
selected_rows = df.iloc[np.linspace(0, len(df) - 1, 5, dtype=int)]

def get_distribution(mean, var, skew):
    std = np.sqrt(var)
    if abs(skew) < 1e-6:
        return norm(loc=mean, scale=std), "Normal"
    else:
        return pearson3(skew=skew, loc=mean, scale=std), "Pearson Type III"

def plot_distributions(rows, direction, output_path):
    plt.figure(figsize=(10, 6))
    colors = plt.cm.viridis(np.linspace(0, 1, len(rows)))

    for idx, row in enumerate(rows.itertuples()):
        if direction == "X":
            mean, var, skew = row.meanX, row.varX, row.skewX
        else:
            mean, var, skew = row.meanY, row.varY, row.skewY

        radius = row.radius
        dist, dist_type = get_distribution(mean, var, skew)

        try:
            x = np.linspace(dist.ppf(0.001), dist.ppf(0.999), 500)
            y = dist.pdf(x)
            plt.plot(x, y, label=f'Radius {radius} ({dist_type})', color=colors[idx])
            plt.axvline(mean, linestyle='--', color=colors[idx], alpha=0.7)
        except Exception as e:
            print(f"Error at radius {radius}: {e}")

    plt.title(f"Approximated Distributions ({direction}-direction)")
    plt.xlabel(direction.lower())
    plt.ylabel("Density")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()

# Generate both plots
plot_distributions(selected_rows, "X", "pearson_distributions_X.png")
plot_distributions(selected_rows, "Y", "pearson_distributions_Y.png")