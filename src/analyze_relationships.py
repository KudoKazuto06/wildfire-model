
from pathlib import Path
import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ============================================================
# 1. CONFIGURATION
# ============================================================

ROOT = Path(__file__).resolve().parent.parent

DATA_PATH = (
    ROOT
    / "data"
    / "raw"
    / "fp-historical-wildfire-data-2006-2025.xlsx"
)

FIGURE_DIR = ROOT / "figures"
FIGURE_DIR.mkdir(parents=True, exist_ok=True)

# Suppress harmless Excel workbook-style warnings
warnings.filterwarnings(
    "ignore",
    message="Workbook contains no default style"
)

# ============================================================
# 2. LOAD DATASET
# ============================================================

print("=" * 65)
print("STAT 471 - WILDFIRE DATA ANALYSIS")
print("=" * 65)

df = pd.read_excel(DATA_PATH)

print(f"\nOriginal observations: {len(df):,}")
print(f"Total variables: {len(df.columns)}")

# ============================================================
# 3. SELECT RELEVANT VARIABLES
# ============================================================

columns = [
    "YEAR",
    "FIRE_NUMBER",
    "FIRE_SPREAD_RATE",
    "WIND_SPEED",
    "RELATIVE_HUMIDITY",
    "TEMPERATURE",
    "FUEL_TYPE",
    "CURRENT_SIZE",
]

data = df[columns].copy()

numeric_columns = [
    "FIRE_SPREAD_RATE",
    "WIND_SPEED",
    "RELATIVE_HUMIDITY",
    "TEMPERATURE",
    "CURRENT_SIZE",
]

for col in numeric_columns:
    data[col] = pd.to_numeric(data[col], errors="coerce")

# Remove rows with missing main numerical variables
data = data.dropna(
    subset=[
        "FIRE_SPREAD_RATE",
        "WIND_SPEED",
        "RELATIVE_HUMIDITY",
        "TEMPERATURE",
    ]
)

# Remove physically invalid observations
data = data[
    (data["FIRE_SPREAD_RATE"] >= 0)
    & (data["WIND_SPEED"] >= 0)
    & (data["RELATIVE_HUMIDITY"] >= 0)
    & (data["RELATIVE_HUMIDITY"] <= 100)
]

print(f"Valid observations: {len(data):,}")
print(f"Removed observations: {len(df) - len(data):,}")

# ============================================================
# 4. DESCRIPTIVE STATISTICS
# ============================================================

print("\n" + "=" * 65)
print("DESCRIPTIVE STATISTICS")
print("=" * 65)

print(
    data[numeric_columns]
    .describe(percentiles=[0.25, 0.5, 0.75, 0.9, 0.95, 0.99])
    .round(3)
    .to_string()
)

print("\nZero spread-rate observations:")
zero_count = (data["FIRE_SPREAD_RATE"] == 0).sum()

print(
    f"{zero_count:,} "
    f"({100 * zero_count / len(data):.2f}%)"
)

# ============================================================
# 5. CORRELATION ANALYSIS
# ============================================================

print("\n" + "=" * 65)
print("PEARSON CORRELATION MATRIX")
print("=" * 65)

correlation_columns = [
    "FIRE_SPREAD_RATE",
    "WIND_SPEED",
    "RELATIVE_HUMIDITY",
    "TEMPERATURE",
]

pearson_corr = data[correlation_columns].corr(
    method="pearson"
)

print(pearson_corr.round(4).to_string())

print("\nSPEARMAN CORRELATION MATRIX")
print("=" * 65)

spearman_corr = data[correlation_columns].corr(
    method="spearman"
)

print(spearman_corr.round(4).to_string())

# ============================================================
# 6. SCATTER PLOTS
# ============================================================

predictors = [
    "WIND_SPEED",
    "RELATIVE_HUMIDITY",
    "TEMPERATURE",
]

for predictor in predictors:

    fig, ax = plt.subplots(figsize=(9, 5))

    ax.scatter(
        data[predictor],
        data["FIRE_SPREAD_RATE"],
        alpha=0.2,
        s=10,
    )

    ax.set_xlabel(predictor.replace("_", " ").title())
    ax.set_ylabel("Recorded Fire Spread Rate")

    ax.set_title(
        f"Fire Spread Rate vs "
        f"{predictor.replace('_', ' ').title()}"
    )

    ax.grid(alpha=0.2)

    fig.tight_layout()

    filename = (
        f"scatter_{predictor.lower()}.png"
    )

    fig.savefig(
        FIGURE_DIR / filename,
        dpi=200
    )

    plt.close(fig)

    print(f"Saved: {filename}")

# ============================================================
# 7. WIND SPEED INTERVAL ANALYSIS
# ============================================================

print("\n" + "=" * 65)
print("FIRE SPREAD RATE BY WIND SPEED INTERVAL")
print("=" * 65)

wind_data = data.copy()

wind_bins = [
    0, 5, 10, 15, 20, 30, 40, 60, np.inf
]

wind_data["WIND_BIN"] = pd.cut(
    wind_data["WIND_SPEED"],
    bins=wind_bins,
    include_lowest=True
)

wind_summary = (
    wind_data
    .groupby("WIND_BIN", observed=True)["FIRE_SPREAD_RATE"]
    .agg(
        count="count",
        mean="mean",
        median="median",
        p90=lambda x: x.quantile(0.90),
        p95=lambda x: x.quantile(0.95),
    )
)

print(wind_summary.round(3).to_string())

# Plot binned statistics

fig, ax = plt.subplots(figsize=(10, 5))

x = np.arange(len(wind_summary))

ax.plot(
    x,
    wind_summary["median"],
    marker="o",
    label="Median"
)

ax.plot(
    x,
    wind_summary["mean"],
    marker="^",
    label="Mean"
)

ax.plot(
    x,
    wind_summary["p90"],
    marker="s",
    label="90th Percentile"
)

ax.set_xticks(x)

ax.set_xticklabels(
    wind_summary.index.astype(str),
    rotation=35
)

ax.set_xlabel("Wind Speed Interval")
ax.set_ylabel("Recorded Fire Spread Rate")
ax.set_title("Fire Spread Rate by Wind Speed Interval")

ax.legend()
ax.grid(alpha=0.2)

fig.tight_layout()

fig.savefig(
    FIGURE_DIR / "wind_speed_binned_analysis.png",
    dpi=200
)

plt.close(fig)

print("\nSaved: wind_speed_binned_analysis.png")

# ============================================================
# 8. FIRE SPREAD RATE DISTRIBUTION
# ============================================================

fig, ax = plt.subplots(figsize=(9, 5))

# log1p handles zero values
log_spread = np.log1p(data["FIRE_SPREAD_RATE"])

ax.hist(
    log_spread,
    bins=50,
    edgecolor="black",
    alpha=0.8
)

ax.set_xlabel("log(1 + Fire Spread Rate)")
ax.set_ylabel("Frequency")

ax.set_title(
    "Distribution of Log-Transformed Fire Spread Rate"
)

ax.grid(alpha=0.2)

fig.tight_layout()

fig.savefig(
    FIGURE_DIR / "fire_spread_distribution.png",
    dpi=200
)

plt.close(fig)

print("Saved: fire_spread_distribution.png")

# ============================================================
# 9. FUEL TYPE ANALYSIS
# ============================================================

print("\n" + "=" * 65)
print("FIRE SPREAD RATE BY FUEL TYPE")
print("=" * 65)

fuel_data = data.dropna(
    subset=["FUEL_TYPE"]
).copy()

fuel_summary = (
    fuel_data
    .groupby("FUEL_TYPE")["FIRE_SPREAD_RATE"]
    .agg(
        count="count",
        mean="mean",
        median="median",
        p90=lambda x: x.quantile(0.90),
    )
    .sort_values("count", ascending=False)
)

print(fuel_summary.round(3).head(15).to_string())

# Plot only the 10 most common fuel types
top_fuels = fuel_summary.head(10).index

plot_data = fuel_data[
    fuel_data["FUEL_TYPE"].isin(top_fuels)
]

fuel_p90 = (
    plot_data
    .groupby("FUEL_TYPE")["FIRE_SPREAD_RATE"]
    .quantile(0.90)
    .sort_values()
)

fig, ax = plt.subplots(figsize=(10, 6))

ax.barh(
    fuel_p90.index.astype(str),
    fuel_p90.values
)

ax.set_xlabel("90th Percentile of Fire Spread Rate")
ax.set_ylabel("Fuel Type")

ax.set_title(
    "High Fire Spread Rates by Fuel Type"
)

ax.grid(axis="x", alpha=0.2)

fig.tight_layout()

fig.savefig(
    FIGURE_DIR / "fuel_type_comparison.png",
    dpi=200
)

plt.close(fig)

print("Saved: fuel_type_comparison.png")

# ============================================================
# 10. SAVE CLEANED DATA
# ============================================================

PROCESSED_DIR = ROOT / "data" / "processed"
PROCESSED_DIR.mkdir(parents=True, exist_ok=True)

output_path = PROCESSED_DIR / "wildfire_cleaned.csv"

data.to_csv(output_path, index=False)

print("\n" + "=" * 65)
print("ANALYSIS COMPLETE")
print("=" * 65)

print(f"Cleaned dataset: {output_path}")
print(f"Figures saved in: {FIGURE_DIR}")
print(f"Total observations analyzed: {len(data):,}")


print("\nFIRE SPREAD RATE FREQUENCIES:")
print(
    data["FIRE_SPREAD_RATE"]
    .value_counts()
    .head(20)
)

print("\nWIND SPEED BIN COUNTS:")
print(wind_summary["count"])

print("\nFUEL TYPE COUNTS:")
print(
    data["FUEL_TYPE"]
    .value_counts()
    .head(15)
)

print("\nPROPORTION OF ZERO SPREAD BY FUEL TYPE:")
print(
    data.groupby("FUEL_TYPE")["FIRE_SPREAD_RATE"]
    .apply(lambda x: (x == 0).mean())
    .sort_values(ascending=False)
    .head(15)
)