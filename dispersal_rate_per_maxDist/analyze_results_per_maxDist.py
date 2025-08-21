import os
import glob
import pandas as pd
import matplotlib.pyplot as plt

def extract_and_save_dispersal_stats(burn_in_value=1000000,
                                     out_summary="beast_data.tsv",
                                     out_grouped="beast_data_grouped.tsv"):
    results = []

    for directory in os.listdir("."):
        if not os.path.isdir(directory):
            continue

        os.chdir(directory)
        try:
            parts = directory.split("_")
            try:
                disp_dist_max = int(parts[0][4:])   # e.g., 'dist10' -> 10
                sim_num = int(parts[1][3:])  # e.g., 'sim3' -> 3
            except (IndexError, ValueError):
                print(f"Skipping directory with unexpected name: {directory}")
                continue

            log_files = glob.glob("*.log")
            if not log_files:
                continue

            df = pd.read_csv(log_files[0], comment='#', delimiter='\t')
            filtered_df = df[df["state"] > burn_in_value]

            # Delete .trees file to save space
            tree_files = glob.glob("*.trees")
            for tf in tree_files:
                os.remove(tf)

            if not filtered_df.empty:
                results.append({
                    "disp_dist_max": disp_dist_max,
                    "simulation": sim_num,
                    "mean": filtered_df["coordinates.diffusionRate"].mean(),
                    "stdev": filtered_df["coordinates.diffusionRate"].std()
                })
        finally:
            os.chdir("..")

    global_df = pd.DataFrame(results)
    global_df.to_csv(out_summary, sep="\t", index=False)

    summary_df = (
        global_df.groupby("disp_dist_max")["mean"]
        .agg([
            ("mean_diffusion_rate", "mean"),
            ("median_diffusion_rate", "median"),
            ("quantile_2.5", lambda x: x.quantile(0.025)),
            ("quantile_97.5", lambda x: x.quantile(0.975))
        ])
        .reset_index()
    )
    summary_df.to_csv(out_grouped, sep="\t", index=False)
    return summary_df

extract_and_save_dispersal_stats()