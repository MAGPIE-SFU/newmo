import os

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages


def create_hist_plot_mass(df):

    print("help me..")

    return None


def sa_observed_histograms(df, overwrite=False):
    output_file = "../plots/sa_observed_histograms.pdf"
    if os.path.exists(output_file) and not overwrite:
        return

    with PdfPages(output_file) as pdf:
        for group, d in df.groupby("group"):
            d = d.sort_values(by="year")

            fig, ax = plt.subplots(figsize=(5, 5))

            sns.barplot(
                data=d,
                x="year",
                y="frac",
                order=sorted(d["year"].unique()),
                ax=ax,
                errorbar=None
            )

            ax.set_title(f"Group: {group}")
            ax.set_ylabel("Fraction")
            ax.set_xlabel("")  # No need to add Year, should be obvious

            ax.tick_params(axis="x", labelrotation=45)

            plt.tight_layout()

            pdf.savefig(fig)
            plt.close(fig)


if __name__ == '__main__':
    # Massachusetts data

    # todo WIP
    # mass_df = pd.read_csv("../data/USA-Mass/")
    # create_hist_plot_mass(mass_df)

    # South Africa data
    sa_real_data = pd.read_csv("../data/SouthAfrica/xt_grouped.csv")
    sa_observed_histograms(sa_real_data)

    # todo model output data

