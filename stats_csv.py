import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from sklearn.decomposition import PCA
from sklearn.impute import SimpleImputer

# Load the data
region1_data = pd.read_csv('Derivation_dataset.csv')
region2_data = pd.read_csv('External_validation_dataset.csv')

# Ensure the columns match between the two datasets
common_features = list(set(region1_data.columns) & set(region2_data.columns))
region1_data = region1_data[common_features]
region2_data = region2_data[common_features]


# 1. Basic Statistical Comparison
def compare_basic_stats(df1, df2, region1_name="Region 1", region2_name="Region 2"):
    results = []

    for feature in df1.columns:
        r1_mean = df1[feature].mean()
        r2_mean = df2[feature].mean()
        r1_var = df1[feature].var()
        r2_var = df2[feature].var()
        r1_median = df1[feature].median()
        r2_median = df2[feature].median()

        # Perform t-test
        t_stat, p_value = stats.ttest_ind(
            df1[feature].dropna(),
            df2[feature].dropna(),
            equal_var=False  # Welch's t-test (doesn't assume equal variances)
        )

        # Perform F-test for variance comparison
        f_stat = r1_var / r2_var if r1_var > r2_var else r2_var / r1_var
        f_pvalue = 1 - stats.f.cdf(f_stat,
                                   df1[feature].dropna().count() - 1,
                                   df2[feature].dropna().count() - 1)

        # KS test for distribution comparison
        ks_stat, ks_pvalue = stats.ks_2samp(
            df1[feature].dropna(),
            df2[feature].dropna()
        )

        results.append({
            'Feature': feature,
            f'{region1_name}_Mean': r1_mean,
            f'{region2_name}_Mean': r2_mean,
            'Mean_Diff': r1_mean - r2_mean,
            'Mean_Diff_Percent': ((r1_mean - r2_mean) / r1_mean * 100) if r1_mean != 0 else np.nan,
            f'{region1_name}_Variance': r1_var,
            f'{region2_name}_Variance': r2_var,
            'Variance_Ratio': r1_var / r2_var if r2_var != 0 else np.nan,
            f'{region1_name}_Median': r1_median,
            f'{region2_name}_Median': r2_median,
            'T_stat': t_stat,
            'T_pvalue': p_value,
            'F_stat': f_stat,
            'F_pvalue': f_pvalue,
            'KS_stat': ks_stat,
            'KS_pvalue': ks_pvalue
        })

    return pd.DataFrame(results)


stats_comparison = compare_basic_stats(region1_data, region2_data)

# Print top features with significant differences
alpha = 0.05
significant_diff_features = stats_comparison[stats_comparison['T_pvalue'] < alpha]
print(f"Features with significantly different means (p < {alpha}):")
print(significant_diff_features[['Feature', 'Region 1_Mean', 'Region 2_Mean', 'T_pvalue']])


# 2. Visualization functions
def create_comparison_visualizations(df1, df2, feature, region1_name="Region 1", region2_name="Region 2"):
    plt.figure(figsize=(15, 5))

    # Histogram
    plt.subplot(1, 3, 1)
    plt.hist(df1[feature].dropna(), alpha=0.5, label=region1_name)
    plt.hist(df2[feature].dropna(), alpha=0.5, label=region2_name)
    plt.title(f'Histogram of {feature}')
    plt.legend()

    # Boxplot (Fix the deprecation warning)
    plt.subplot(1, 3, 2)
    boxplot_data = [df1[feature].dropna(), df2[feature].dropna()]
    plt.boxplot(boxplot_data, tick_labels=[region1_name, region2_name])  # Changed labels to tick_labels
    plt.title(f'Boxplot of {feature}')

    # QQ Plot
    plt.subplot(1, 3, 3)
    stats.probplot(df1[feature].dropna(), plot=plt)
    plt.title(f'QQ Plot of {feature} for {region1_name}')

    plt.tight_layout()
    plt.savefig(f'comparison_{feature}.png')
    plt.close()


# Create visualizations for top 5 most different features
for feature in significant_diff_features['Feature'].head(5):
    create_comparison_visualizations(region1_data, region2_data, feature)


# 3. Correlation Matrix Comparison
def compare_correlation_matrices(df1, df2, region1_name="Region 1", region2_name="Region 2"):
    corr1 = df1.corr()
    corr2 = df2.corr()

    # Calculate the difference between correlation matrices
    corr_diff = corr1 - corr2

    # Visualize the correlation matrices and their difference
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))

    sns.heatmap(corr1, ax=axes[0], cmap='coolwarm', vmin=-1, vmax=1)
    axes[0].set_title(f'Correlation Matrix - {region1_name}')

    sns.heatmap(corr2, ax=axes[1], cmap='coolwarm', vmin=-1, vmax=1)
    axes[1].set_title(f'Correlation Matrix - {region2_name}')

    sns.heatmap(corr_diff, ax=axes[2], cmap='coolwarm', vmin=-1, vmax=1)
    axes[2].set_title('Correlation Difference (Region1 - Region2)')

    plt.tight_layout()
    plt.savefig('correlation_comparison.png')
    plt.close()

    return corr1, corr2, corr_diff


corr1, corr2, corr_diff = compare_correlation_matrices(region1_data, region2_data)

# Find feature pairs with the largest correlation differences
corr_diff_df = corr_diff.unstack().reset_index()
corr_diff_df.columns = ['Feature1', 'Feature2', 'Correlation_Diff']
corr_diff_df = corr_diff_df[corr_diff_df['Feature1'] != corr_diff_df['Feature2']]
corr_diff_df = corr_diff_df.sort_values(by='Correlation_Diff', key=abs, ascending=False)
print("Feature pairs with largest correlation differences:")
print(corr_diff_df.head(10))


# 4. PCA for dimensionality reduction and visualization (Fixed to handle NaN values)
def perform_pca_comparison(df1, df2, region1_name="Region 1", region2_name="Region 2"):
    # Combine datasets for PCA
    df1_copy = df1.copy()
    df2_copy = df2.copy()

    df1_copy['Region'] = region1_name
    df2_copy['Region'] = region2_name

    combined_data = pd.concat([df1_copy, df2_copy], ignore_index=True)

    # Extract features and region labels
    features = [col for col in combined_data.columns if col != 'Region']
    X = combined_data[features]
    regions = combined_data['Region']

    # Handle missing values using mean imputation
    imputer = SimpleImputer(strategy='mean')
    X_imputed = imputer.fit_transform(X)

    # Standardize the data
    X_std = (X_imputed - np.mean(X_imputed, axis=0)) / np.std(X_imputed, axis=0)

    # Apply PCA
    pca = PCA(n_components=2)  # Reduce to 2 dimensions for visualization
    X_pca = pca.fit_transform(X_std)

    # Create PCA dataframe for plotting
    pca_df = pd.DataFrame({
        'PC1': X_pca[:, 0],
        'PC2': X_pca[:, 1],
        'Region': regions.values
    })

    # Plot PCA results
    plt.figure(figsize=(10, 8))

    for region, color in zip([region1_name, region2_name], ['blue', 'red']):
        subset = pca_df[pca_df['Region'] == region]
        plt.scatter(subset['PC1'], subset['PC2'], c=color, label=region, alpha=0.6)

    plt.xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.2%} variance)')
    plt.ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.2%} variance)')
    plt.title('PCA of ECG Features by Region')
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.7)

    plt.savefig('pca_comparison.png')
    plt.close()

    return pca, pca_df


pca, pca_df = perform_pca_comparison(region1_data, region2_data)


# 5. Generate a comprehensive report
def generate_report(stats_df, pca_var_ratio):
    print("\n===== ECG FEATURE COMPARISON REPORT =====\n")

    # Overall similarity assessment
    significant_features = stats_df[stats_df['T_pvalue'] < 0.05]
    pct_significant = len(significant_features) / len(stats_df) * 100

    print(
        f"OVERALL SIMILARITY: {100 - pct_significant:.1f}% of features show no significant difference between regions")

    # Most different features
    print("\nMOST DIFFERENT FEATURES (by t-test p-value):")
    most_diff = stats_df.sort_values('T_pvalue').head(5)
    for idx, row in most_diff.iterrows():
        print(
            f"- {row['Feature']}: Region 1 mean = {row['Region 1_Mean']:.3f}, Region 2 mean = {row['Region 2_Mean']:.3f}, p-value = {row['T_pvalue']:.5f}")

    # Most similar features
    print("\nMOST SIMILAR FEATURES (by t-test p-value):")
    most_similar = stats_df.sort_values('T_pvalue', ascending=False).head(5)
    for idx, row in most_similar.iterrows():
        print(
            f"- {row['Feature']}: Region 1 mean = {row['Region 1_Mean']:.3f}, Region 2 mean = {row['Region 2_Mean']:.3f}, p-value = {row['T_pvalue']:.5f}")

    # PCA information
    print(f"\nPCA ANALYSIS: First two components explain {(pca_var_ratio[0] + pca_var_ratio[1]):.2%} of total variance")

    print("\n========================================\n")


generate_report(stats_comparison, pca.explained_variance_ratio_)


# 6. Additional analysis: Missing data comparison
def compare_missing_data(df1, df2, region1_name="Region 1", region2_name="Region 2"):
    missing_stats = []

    for feature in df1.columns:
        r1_missing = df1[feature].isna().sum()
        r1_missing_pct = r1_missing / len(df1) * 100

        r2_missing = df2[feature].isna().sum()
        r2_missing_pct = r2_missing / len(df2) * 100

        missing_stats.append({
            'Feature': feature,
            f'{region1_name}_Missing': r1_missing,
            f'{region1_name}_Missing_Pct': r1_missing_pct,
            f'{region2_name}_Missing': r2_missing,
            f'{region2_name}_Missing_Pct': r2_missing_pct,
            'Missing_Diff_Pct': abs(r1_missing_pct - r2_missing_pct)
        })

    missing_df = pd.DataFrame(missing_stats)

    # Print features with highest missing data difference
    print("\nFeatures with highest missing data percentage difference between regions:")
    print(missing_df.sort_values('Missing_Diff_Pct', ascending=False).head(10))

    # Plot missing data comparison for top features
    top_missing_features = missing_df.sort_values('Missing_Diff_Pct', ascending=False).head(5)['Feature'].tolist()

    plt.figure(figsize=(12, 6))
    x = np.arange(len(top_missing_features))
    width = 0.35

    r1_missing = [missing_df[missing_df['Feature'] == f][f'{region1_name}_Missing_Pct'].values[0] for f in
                  top_missing_features]
    r2_missing = [missing_df[missing_df['Feature'] == f][f'{region2_name}_Missing_Pct'].values[0] for f in
                  top_missing_features]

    plt.bar(x - width / 2, r1_missing, width, label=region1_name)
    plt.bar(x + width / 2, r2_missing, width, label=region2_name)

    plt.xlabel('Features')
    plt.ylabel('Missing Data (%)')
    plt.title('Missing Data Comparison Between Regions')
    plt.xticks(x, top_missing_features, rotation=45)
    plt.legend()
    plt.tight_layout()
    plt.savefig('missing_data_comparison.png')
    plt.close()

    return missing_df


missing_comparison = compare_missing_data(region1_data, region2_data)