import os

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
from sklearn.metrics import confusion_matrix, balanced_accuracy_score
from histogram import FIXED_POPULATION_SIZE


def get_trends(values, threshold=0.07):
    trends = []

    scale = max(values) - min(values)

    for prev, curr in zip(values[:-1], values[1:]):
        if scale == 0:
            trends.append("stable")
            continue

        normalized_change = (curr - prev) / scale

        if normalized_change > threshold:
            trends.append("up")
        elif normalized_change < -threshold:
            trends.append("down")
        else:
            trends.append("stable")

    return trends


def compare_sim_obs(observed_df, simulated_df, overwrite=False):
    output_file = "../plots/sa_confusion_matrix_trends.pdf"
    if os.path.exists(output_file) and not overwrite:
        return

    # the value column is counts
    sim_long = simulated_df.melt(
        id_vars="id",
        var_name="month",
        value_name="value"
    )
    sim_long["month"] = sim_long["month"].str.replace("month\\.", "", regex=True)

    sim_long["frac"] = sim_long["value"] / FIXED_POPULATION_SIZE

    all_obs = []
    all_sim = []

    # compute a confusion matrix obsv vs simulated: trends are up down or stable
    for group, obs_d in sorted(observed_df.groupby("group"), key=lambda x: x[0]):
        sim_d = sim_long[sim_long["id"] == group]

        obs_d = obs_d.sort_values("year")
        sim_d = sim_d.sort_values("month")

        obs = obs_d["frac"].to_numpy()
        sim = sim_d["frac"].to_numpy()

        obs_trends = get_trends(obs)
        sim_trends = get_trends(sim)

        all_obs.extend(obs_trends)
        all_sim.extend(sim_trends)

    labels = ["up", "stable", "down"]

    cm = confusion_matrix(
        all_obs,
        all_sim,
        labels=labels
    )

    cm_df = pd.DataFrame(
        cm,
        index=["Observed Up", "Observed Stable", "Observed Down"],
        columns=["Predicted Up", "Predicted Stable", "Predicted Down"]
    )

    plt.figure(figsize=(6, 5))

    sns.heatmap(
        cm_df,
        annot=True,
        fmt="d",
        cmap="Blues",
        cbar=False,
        linewidths=0.5
    )

    plt.title("Trend Prediction Confusion Matrix")
    # plt.ylabel("Observed")
    # plt.xlabel("Simulated")

    plt.tight_layout()
    plt.savefig(output_file, bbox_inches="tight")
    plt.close()


def explore_confusion_matrices(obs_df, simulated_df, overwrite=False):
    output_file = "../plots/exploration_cm.pdf"
    if os.path.exists(output_file) and not overwrite:
        return

    labels = ["up", "stable", "down"]

    # simulated_df is already long format:
    # serotype, month, value, simulation_id, v, sigma, migr
    simulated_df = simulated_df.copy()
    simulated_df["frac"] = simulated_df["value"] / FIXED_POPULATION_SIZE
    # simulated_df["frac"] = simulated_df.groupby(["simulation_id", "month"])["value"].transform(
    #     lambda x: x / x.sum()
    # )

    parameter_cols = ["v", "sigma", "migr"]

    skipped_groups = []

    scores = []  # collect one row per simulation

    with PdfPages(output_file) as pdf:
        for sim_id, sim_group in sorted(
                simulated_df.groupby("simulation_id"), key=lambda x: x[0]
        ):
            params = sim_group[parameter_cols].iloc[0]

            all_obs = []
            all_sim = []

            for group, obs_d in sorted(obs_df.groupby("group"), key=lambda x: x[0]):
                sim_d = sim_group[sim_group["serotype"] == group]
                obs_d = obs_d.sort_values("year")

                if not (obs_d["frac"] >= 0.025).any():
                    skipped_groups.append(group)
                    continue

                sim_d = sim_d.sort_values("month")
                obs = obs_d["frac"].to_numpy()
                sim = sim_d["frac"].to_numpy()

                obs_trends = get_trends(obs)
                sim_trends = get_trends(sim)

                all_obs.extend(obs_trends)
                all_sim.extend(sim_trends)

            cm = confusion_matrix(all_obs, all_sim, labels=labels)

            # --- diagonal-based scoring ---
            n_total = cm.sum()
            n_correct = np.trace(cm)
            accuracy = n_correct / n_total if n_total > 0 else np.nan
            bal_acc = balanced_accuracy_score(all_obs, all_sim) if n_total > 0 else np.nan

            scores.append({
                "simulation_id": sim_id,
                "accuracy": accuracy,
                "balanced_accuracy": bal_acc,
                "n_points": n_total,
                **params.to_dict(),
            })

            cm_df = pd.DataFrame(
                cm,
                index=["Observed Up", "Observed Stable", "Observed Down"],
                columns=["Predicted Up", "Predicted Stable", "Predicted Down"]
            )

            plt.figure(figsize=(6, 5))
            sns.heatmap(cm_df, annot=True, fmt="d", cmap="Blues", cbar=False, linewidths=0.5)
            plt.title(
                f"Simulation {sim_id}: v={params['v']}, sigma={params['sigma']}, "
                f"migr={params['migr']}\nAccuracy={accuracy:.2f}, Balanced={bal_acc:.2f}"
            )
            plt.tight_layout()
            pdf.savefig(bbox_inches="tight")
            plt.close()

    # --- ranking summary, separate PDF ---
    scores_df = pd.DataFrame(scores).sort_values("balanced_accuracy", ascending=False)

    ranking_file = output_file.replace(".pdf", "_ranking.pdf")

    with PdfPages(ranking_file) as pdf:
        # Page 1: bar chart of top N by balanced accuracy
        top_n = min(30, len(scores_df))
        plt.figure(figsize=(8, max(4, top_n * 0.3)))
        plot_df = scores_df.head(top_n).iloc[::-1]  # reverse so best is on top
        plt.barh(plot_df["simulation_id"].astype(str), plot_df["balanced_accuracy"],
                 color="steelblue")
        plt.xlabel("Balanced Accuracy")
        plt.ylabel("Simulation ID")
        plt.title(f"Top {top_n} Simulations by Balanced Accuracy")
        plt.tight_layout()
        pdf.savefig(bbox_inches="tight")
        plt.close()

        # Page 2: full table
        fig, ax = plt.subplots(figsize=(8.5, max(4, len(scores_df) * 0.25)))
        ax.axis("off")
        table_cols = ["simulation_id", "accuracy", "balanced_accuracy", "n_points"] + parameter_cols
        tbl = ax.table(
            cellText=scores_df[table_cols].round(3).values,
            colLabels=table_cols,
            loc="center",
            cellLoc="center",
        )
        tbl.auto_set_font_size(False)
        tbl.set_fontsize(7)
        tbl.scale(1, 1.2)
        plt.title("Full Ranking", pad=20)
        pdf.savefig(bbox_inches="tight")
        plt.close()

    # also handy to have as CSV for later analysis (e.g. feeding into SALib)
    scores_df.to_csv(output_file.replace(".pdf", "_ranking.csv"), index=False)

    # print(f"Skipped groups: {set(skipped_groups)}")


if __name__ == '__main__':
    # South Africa data
    sa_real_data = pd.read_csv("../data/SouthAfrica/xt_grouped.csv")

    # model only starts in 2008, so we delete earlier times
    sa_real_data = sa_real_data[sa_real_data["year"] > 2007]

    sa_simulated_data_exploration = pd.read_csv("../data/SouthAfrica/parameter_exploration.csv")
    sa_simulated_data = pd.read_csv("../data/SouthAfrica/simulation_data.csv")

    # compare_sim_obs(sa_real_data, sa_simulated_data, overwrite=True)
    explore_confusion_matrices(sa_real_data, sa_simulated_data_exploration, overwrite=True)
