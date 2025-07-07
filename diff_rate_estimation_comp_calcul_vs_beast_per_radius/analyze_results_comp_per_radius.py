import os
import glob
import pandas as pd
import matplotlib.pyplot as plt

def extract_empirical_dispersal(full_output="empirical_dispersal_all.tsv",
                                grouped_output="empirical_dispersal_grouped.tsv",
                                moment_output="empirical_dispersal_moments_all.tsv",
                                moment_grouped_output="empirical_dispersal_moments_grouped.tsv"):

    empirical_dispersal_data = []
    empirical_moments_data = []

    for directory in os.listdir("."):
        if not os.path.isdir(directory):
            continue

        os.chdir(directory)
        try:
            parts = directory.split("_")
            radius = int(parts[0][1:])   # e.g., 'r10' -> 10
            sim_num = int(parts[1][3:])  # e.g., 'sim3' -> 3

            for file in os.listdir("."):
                if file.endswith("_GSpace_param_summary.txt"):
                    with open(file, "r") as f:
                        lines = f.readlines()

                    # Extract dispersal distribution
                    start_idx = None
                    for i, line in enumerate(lines):
                        if line.strip().startswith("Empirical dispersal distribution"):
                            start_idx = i + 1
                            break

                    if start_idx is not None:
                        for j in range(start_idx, len(lines)):
                            line = lines[j].strip()
                            if not line or line.startswith("#") or line.startswith("%%"):
                                break
                            parts = line.split()
                            if len(parts) != 3:
                                continue
                            step, freqX, freqY = parts
                            empirical_dispersal_data.append({
                                "radius": radius,
                                "simulation": sim_num,
                                "step": int(step),
                                "freqX": float(freqX),
                                "freqY": float(freqY)
                            })

                    # Extract dispersal moments
                    meanX = meanY = varX = varY = skewX = skewY = kurtX = kurtY = None
                    for line in lines:
                        if line.startswith("Empirical dispersal mean"):
                            meanX, meanY = map(float, line.split(":")[1].strip().split())
                        elif line.startswith("Empirical dispersal sigma2"):
                            varX, varY = map(float, line.split(":")[1].strip().split())
                        elif line.startswith("Empirical dispersal skewness"):
                            skewX, skewY = map(float, line.split(":")[1].strip().split())
                        elif line.startswith("Empirical dispersal kurtosis"):
                            kurtX, kurtY = map(float, line.split(":")[1].strip().split())

                    empirical_moments_data.append({
                        "radius": radius,
                        "simulation": sim_num,
                        "meanX": meanX,
                        "meanY": meanY,
                        "varX": varX,
                        "varY": varY,
                        "skewX": skewX,
                        "skewY": skewY,
                        "kurtX": kurtX,
                        "kurtY": kurtY
                    })
        finally:
            os.chdir("..")

    empirical_df = pd.DataFrame(empirical_dispersal_data)
    empirical_df.to_csv(full_output, sep="\t", index=False)

    empirical_summary = (
        empirical_df.groupby(["radius", "step"])[["freqX", "freqY"]]
        .mean()
        .reset_index()
    )
    empirical_summary.to_csv(grouped_output, sep="\t", index=False)

    moments_df = pd.DataFrame(empirical_moments_data)
    moments_df.to_csv(moment_output, sep="\t", index=False)

    grouped_moments = (
        moments_df.groupby("radius")[["meanX", "meanY", "varX", "varY", "skewX", "skewY", "kurtX", "kurtY"]]
        .mean()
        .reset_index()
    )
    grouped_moments.to_csv(moment_grouped_output, sep="\t", index=False)

    # Plot each moment parameter (mean, var, skew, kurt) in a separate subplot with X and Y
    fig, axs = plt.subplots(2, 2, figsize=(14, 10))
    axs = axs.flatten()
    parameters = [("meanX", "meanY", "Mean"), ("varX", "varY", "Variance"),
                  ("skewX", "skewY", "Skewness"), ("kurtX", "kurtY", "Kurtosis")]

    for i, (x_param, y_param, title) in enumerate(parameters):
        axs[i].plot(grouped_moments["radius"], grouped_moments[x_param], label="X", color="blue")
        axs[i].plot(grouped_moments["radius"], grouped_moments[y_param], label="Y", color="orange")
        axs[i].set_title(f"{title} by Radius")
        axs[i].set_xlabel("Radius")
        axs[i].set_ylabel(title)
        axs[i].legend()
        axs[i].grid(True)

    plt.tight_layout()
    plt.savefig("empirical_dispersal_moments_by_radius.png",dpi=300)

#extract_empirical_dispersal()

def extract_and_compare_diffusion_estimates(
    burn_in_fraction=0.1,
    out_raw="diff_rate_comparison.tsv",
    out_grouped="diff_rate_comparison_grouped.tsv"):
    from io import StringIO

    records = []

    for directory in os.listdir("."):
        if not os.path.isdir(directory):
            continue

        os.chdir(directory)
        try:
            parts = directory.split("_")
            radius = int(parts[0][1:])   # e.g., 'r10' -> 10
            sim_num = int(parts[1][3:])  # e.g., 'sim3' -> 3

            # Load GSpace summary
            gspace_file = glob.glob("*_seq_char.txt")
            if not gspace_file:
                continue
            gspace_df = pd.read_csv(gspace_file[0], sep=r'\s+', header=0)
            var_x = gspace_df["X"].var()

            # Load BEAST log
            log_files = glob.glob("*.log")
            if not log_files:
                continue
            with open(log_files[0]) as f:
                lines = [line for line in f if not line.startswith("Trace") and not line.startswith("#")]
            log_df = pd.read_csv(StringIO("".join(lines)), sep="\t")
            log_df["age(root)"] = pd.to_numeric(log_df["age(root)"], errors="coerce")
            log_df["coordinates.diffusionRate"] = pd.to_numeric(log_df["coordinates.diffusionRate"], errors="coerce")
            burnin = int(burn_in_fraction * len(log_df))
            log_df_post = log_df.iloc[burnin:]

            if not log_df_post.empty:
                mean_age_root = log_df_post["age(root)"].mean()
                mean_diff_rate = log_df_post["coordinates.diffusionRate"].mean()
                calc_diff_rate = var_x / mean_age_root if mean_age_root > 0 else None

                records.append({
                    "radius": radius,
                    "simulation": sim_num,
                    "varX": var_x,
                    "mean_age_root": mean_age_root,
                    "mean_diff_rate_beast": mean_diff_rate,
                    "calc_diff_rate": calc_diff_rate
                })
        finally:
            os.chdir("..")

    df = pd.DataFrame(records)
    df.to_csv(out_raw, sep="\t", index=False)

    grouped = (
        df.groupby("radius")[["mean_diff_rate_beast", "calc_diff_rate"]]
        .agg([
            ("mean", "mean"),
            ("2.5%", lambda x: x.quantile(0.025)),
            ("97.5%", lambda x: x.quantile(0.975))
        ])
    )
    grouped.columns = ['_'.join(col).strip() for col in grouped.columns]
    grouped = grouped.reset_index()
    grouped.to_csv(out_grouped, sep="\t", index=False)

    # Plot comparison
    plt.figure(figsize=(10, 6))
    plt.plot(grouped["radius"], grouped["mean_calc_diff_rate"], marker='o', label="Var(X)/T Estimate")
    plt.plot(grouped["radius"], grouped["mean_mean_diff_rate_beast"], marker='s', label="BEAST Estimate")
    plt.fill_between(
        grouped["radius"],
        grouped["2.5%_calc_diff_rate"],
        grouped["97.5%_calc_diff_rate"],
        color='blue',
        alpha=0.2,
        label="Var(X)/T 95% CI"
    )
    plt.fill_between(
        grouped["radius"],
        grouped["2.5%_mean_diff_rate_beast"],
        grouped["97.5%_mean_diff_rate_beast"],
        color='orange',
        alpha=0.2,
        label="BEAST 95% CI"
    )
    plt.xlabel("Radius")
    plt.ylabel("Diffusion Rate")
    plt.title("Comparison of Diffusion Rate Estimates")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig("diffusion_rate_comparison.png", dpi=300)

extract_and_compare_diffusion_estimates()