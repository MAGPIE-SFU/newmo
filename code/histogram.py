import os

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
from sklearn.metrics import confusion_matrix


FIXED_POPULATION_SIZE = 100000


def create_hist_plot_mass(df):

    print("help me..")

    return None


def _bar_plot_page(df, x, y, title, ax, x_order=None):
    sns.barplot(
        data=df,
        x=x,
        y=y,
        ax=ax,
        order=x_order,
        errorbar=None
    )

    ax.set_title(title)
    ax.set_ylabel("Fraction")
    ax.set_xlabel("")
    ax.tick_params(axis="x", labelrotation=45)


def _write_grouped_pdf(df, group_col, output_file, x, y, title_prefix, x_order=None, overwrite=False):
    if os.path.exists(output_file) and not overwrite:
        return

    with PdfPages(output_file) as pdf:
        for group, d in df.groupby(group_col):
            d = d.copy()

            fig, ax = plt.subplots(figsize=(5, 5))

            _bar_plot_page(
                d,
                x=x,
                y=y,
                title=f"{title_prefix}: {group}",
                ax=ax,
                x_order=x_order
            )

            plt.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)


def sa_combined_pdf(observed_df, simulated_df, overwrite=False):
    output_file = "../plots/sa_combined_histograms.pdf"
    if os.path.exists(output_file) and not overwrite:
        return

    sim_long = simulated_df.melt(
        id_vars="id",
        var_name="month",
        value_name="value"
    )
    sim_long["month"] = sim_long["month"].str.replace("month\\.", "", regex=True)

    sim_long["frac"] = sim_long["value"] / FIXED_POPULATION_SIZE

    month_order = sorted(sim_long["month"].unique(), key=int)

    with PdfPages(output_file) as pdf:

        for group, obs_d in sorted(observed_df.groupby("group"), key=lambda x: x[0]):

            sim_d = sim_long[sim_long["id"] == group]

            fig, (ax1, ax2) = plt.subplots(
                1, 2,
                figsize=(10, 4),
                sharey=False
            )

            fig.suptitle(f"Group: {group}")

            # -----------------------
            # Observed plot
            # -----------------------
            sns.barplot(
                data=obs_d,
                x="year",
                y="frac",
                order=sorted(obs_d["year"].unique()),
                ax=ax1,
                errorbar=None
            )
            ax1.set_title(f"Observed")
            ax1.set_xlabel("Year")
            ax1.set_ylabel("Fraction")
            ax1.tick_params(axis="x", rotation=45)

            # -----------------------
            # Simulated plot
            # -----------------------
            sns.barplot(
                data=sim_d,
                x="month",
                y="frac",
                order=month_order,
                ax=ax2,
                errorbar=None
            )
            ax2.set_title(f"Simulated")
            ax2.set_xlabel("Month")
            ax2.set_ylabel("Fraction")
            ax2.tick_params(axis="x", rotation=45)

            plt.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)


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


def get_trends(values, threshold=0.05):
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


def population_size_model_plot(sim_df):
    sim_long = sim_df.melt(
        id_vars="id",
        var_name="month",
        value_name="value"
    )
    sim_long["month"] = sim_long["month"].str.replace("month\\.", "", regex=True)

    monthly_population_simulation = sim_long.groupby("month")["value"].sum().reset_index()

    fig, ax = plt.subplots(figsize=(10, 4))

    ax.plot(monthly_population_simulation.month, monthly_population_simulation.value, marker="o")
    ax.axhline(y=FIXED_POPULATION_SIZE, color="red", linestyle="dashed")

    ax.set_xlabel("Month")
    ax.set_ylabel("Population")

    plt.tight_layout()
    plt.savefig("../plots/population_size_model.pdf")
    plt.close()


if __name__ == '__main__':
    # Massachusetts data

    # todo WIP
    # mass_df = pd.read_csv("../data/USA-Mass/")
    # create_hist_plot_mass(mass_df)

    # South Africa data
    sa_real_data = pd.read_csv("../data/SouthAfrica/xt_grouped.csv")

    # model only starts in 2008, so we delete earlier times
    sa_real_data = sa_real_data[sa_real_data["year"] > 2007]

    sa_simulated_data = pd.read_csv("../data/SouthAfrica/simulation_data.csv")

    # one pdf with both bar plots
    # sa_combined_pdf(sa_real_data, sa_simulated_data)
    # compare_sim_obs(sa_real_data, sa_simulated_data)

    population_size_model_plot(sa_simulated_data)
