import os

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages


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


def sa_observed_histograms(df, overwrite=False):
    output_file = "../plots/sa_observed_histograms.pdf"

    _write_grouped_pdf(
        df=df,
        group_col="group",
        output_file=output_file,
        x="year",
        y="frac",
        title_prefix="Group",
        overwrite=overwrite
    )


def sa_simulated_histograms(df, overwrite=False):
    output_file = "../plots/sa_simulated_histograms.pdf"

    if os.path.exists(output_file) and not overwrite:
        return

    df_long = df.melt(
        id_vars="id",
        var_name="month",
        value_name="value"
    )

    df_long["month"] = df_long["month"].str.replace("month\\.", "", regex=True)

    month_order = sorted(df_long["month"].unique(), key=int)

    _write_grouped_pdf(
        df=df_long,
        group_col="id",
        output_file=output_file,
        x="month",
        y="value",
        title_prefix="ID",
        x_order=month_order,
        overwrite=True
    )


if __name__ == '__main__':
    # Massachusetts data

    # todo WIP
    # mass_df = pd.read_csv("../data/USA-Mass/")
    # create_hist_plot_mass(mass_df)

    # South Africa data
    sa_real_data = pd.read_csv("../data/SouthAfrica/xt_grouped.csv")
    sa_observed_histograms(sa_real_data)

    sa_simulated_data = pd.read_csv("../data/SouthAfrica/simulation_data.csv")
    sa_simulated_histograms(sa_simulated_data)
