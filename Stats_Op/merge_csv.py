import pandas as pd
import os


def merge_csvs_same_columns(csv1_path, csv2_path, output_path):
    """
    Append the second CSV file to the first one, column by column.
    Each column in the first CSV will have the corresponding column from the second CSV
    appended to it (excluding the header row of the second CSV).
    """
    # Read the first CSV
    df1 = pd.read_csv(csv1_path)

    # Read the second CSV, skipping the first row
    df2 = pd.read_csv(csv2_path, skiprows=1, header=None)

    # Make sure both have the same number of columns
    if len(df1.columns) != len(df2.columns):
        min_cols = min(len(df1.columns), len(df2.columns))
        df1 = df1.iloc[:, :min_cols]
        df2 = df2.iloc[:, :min_cols]
        print(f"Warning: Using only the first {min_cols} columns from both files")

    # Append each column of df2 to the corresponding column of df1
    result_df = pd.DataFrame()

    for i, col_name in enumerate(df1.columns):
        # Combine the values from both dataframes for this column
        combined_values = pd.concat([df1.iloc[:, i], df2.iloc[:, i]]).reset_index(drop=True)
        result_df[col_name] = combined_values

    # Write the merged dataframe to a new CSV
    result_df.to_csv(output_path, index=False)

    print(f"Successfully merged CSVs with appended columns into {output_path}")
    return result_df


merged_data = merge_csvs_same_columns("Derivation_dataset.csv", "External_validation_dataset.csv", "merged_file.csv")