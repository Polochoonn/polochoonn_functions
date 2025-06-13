import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from typing import Optional, Union

def plot_forest(
    df: Optional[Union[str, pd.DataFrame]] = None,
    csv_path: Optional[str] = None,
    output_path: str = "forest_plot.png",
    top_n: int = 20,
    odds_ratio_col: str = "odds_ratio",
    ci_low_col: str = "or_ci_low",
    ci_high_col: str = "or_ci_high",
    label_col: str = "metabolite",
    pval_col: str = "pval",
    fdr_col: str = "pval_adj",
    plot_title: str = "Top Associations",
    xlabel: str = "Odds Ratio (95% CI)",
    figsize: tuple = (15, 8),
    show_pvals: bool = True,
    fontsize: int = 11
):
    """
    Create a forest plot from a regression results DataFrame or CSV.

    Parameters
    ----------
    df : pd.DataFrame or None
        DataFrame of regression results. If None, uses csv_path.
    csv_path : str or None
        Path to CSV if df is not provided.
    output_path : str
        Filepath for saving the plot (e.g., 'output.png').
    top_n : int
        Number of top results to plot.
    odds_ratio_col, ci_low_col, ci_high_col, label_col, pval_col, fdr_col : str
        Column names for odds ratio, CIs, label, p-value, FDR.
    plot_title : str
        Title for the plot.
    xlabel : str
        X-axis label.
    figsize : tuple
        Figure size (width, height).
    show_pvals : bool
        Whether to annotate p-value/FDR.
    fontsize : int
        Y-label font size.
    """
    # Load data if only CSV path provided
    if df is None and csv_path is not None:
        df = pd.read_csv(csv_path)
    if df is None:
        raise ValueError("Provide either a DataFrame (df) or csv_path.")

    # Sort by FDR or p-value
    if fdr_col in df.columns:
        df = df.sort_values(fdr_col)
        pval_to_show = fdr_col
    else:
        df = df.sort_values(pval_col)
        pval_to_show = pval_col

    # Select top N and reverse for forest plot style
    df_top = df.head(top_n).copy().iloc[::-1]

    # Compute error bars for odds ratio
    error_lower = df_top[odds_ratio_col] - df_top[ci_low_col]
    error_upper = df_top[ci_high_col] - df_top[odds_ratio_col]

    fig, ax = plt.subplots(figsize=figsize)

    # Plot each point with error bars and annotation
    for i, (odds, err_l, err_u, name, ci_l, ci_h, p) in enumerate(
        zip(
            df_top[odds_ratio_col],
            error_lower,
            error_upper,
            df_top[label_col],
            df_top[ci_low_col],
            df_top[ci_high_col],
            df_top[pval_to_show]
        )
    ):
        ax.errorbar(
            odds, i, xerr=[[err_l], [err_u]],
            fmt='o', color='black', capsize=4
        )
        annot = f"OR={odds:.2f}\n({ci_l:.2f}-{ci_h:.2f})"
        if show_pvals:
            if pval_to_show == fdr_col:
                annot += f"\nFDR={p:.3g}"
            else:
                annot += f"\np={p:.3g}"
        ax.text(odds + err_u + 0.05, i, annot, va='center', fontsize=fontsize-1, color='black')

    ax.set_yticks(range(len(df_top)))
    ax.set_yticklabels(df_top[label_col], fontsize=fontsize)
    ax.axvline(1, color='grey', linestyle='--', linewidth=1)
    ax.set_xlabel(xlabel, fontsize=fontsize)
    ax.set_title(plot_title, fontsize=fontsize+2)
    plt.gca().invert_yaxis()
    plt.tight_layout()
    plt.savefig(output_path, dpi=200)
    plt.close()
    print(f"Forest plot saved to {output_path}")

# Example usage:
# plot_forest(csv_path="path/to/logistic_regression_ALL.csv", output_path="top_20_metabolite_associations.png", top_n=20, plot_title="Top 20 Metabolite Associations with CV Death/MI")