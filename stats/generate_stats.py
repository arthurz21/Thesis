import pandas as pd
import numpy as np


def add_statistics_to_csv(input_file, output_file):
    """
    Read CSV file, calculate statistics for each numeric column,
    and append results as a new row.

    Args:
        input_file (str): Path to input CSV file
        output_file (str): Path to save output CSV file
    """
    # Read the CSV file
    df = pd.read_csv(input_file)

    # Initialize dictionaries for statistics
    means = {}
    stds = {}

    # Calculate statistics for numeric columns
    for column in df.columns:
        if pd.api.types.is_numeric_dtype(df[column]):
            means[column] = df[column].mean()
            stds[column] = df[column].std()
        else:
            means[column] = ''
            stds[column] = ''

    # Create new rows for statistics
    mean_row = pd.DataFrame([means], index=['Average'])
    std_row = pd.DataFrame([stds], index=['Std Dev'])

    # Concatenate original data with statistics rows
    result = pd.concat([df, mean_row, std_row])

    # Save to new CSV file
    result.to_csv(output_file, index=True)

# Example usage
add_statistics_to_csv('Derivation_dataset.csv', 'output_with_stats.csv')