import pandas as pd
import numpy as np


def add_numeric_statistics_to_csv(input_file, output_file):
    """
    Read CSV file, calculate statistics for numeric columns only,
    and append results as new rows. Text columns are preserved but ignored for statistics.
    Results are rounded to 6 decimal places.

    Args:
        input_file (str): Path to input CSV file
        output_file (str): Path to save output CSV file
    """
    # Read the CSV file
    df = pd.read_csv(input_file)

    # Initialize dictionaries for statistics
    means = {}
    stds = {}

    # Calculate statistics only for numeric columns
    for column in df.columns:
        # Check if column contains numeric data
        if df[column].dtype in ['int64', 'float64']:
            means[column] = round(df[column].mean(), 6)
            stds[column] = round(df[column].std(), 6)
        else:
            # For text columns, use empty string
            means[column] = ''
            stds[column] = ''

    # Create new rows for statistics
    mean_row = pd.DataFrame([means], index=['Average'])
    std_row = pd.DataFrame([stds], index=['Std Dev'])

    # Concatenate original data with statistics rows
    result = pd.concat([df, mean_row, std_row])

    # Save to new CSV file
    result.to_csv(output_file, index=True)

    # Print summary of which columns were processed
    print("\nStatistics calculated for these numeric columns:")
    for column in df.columns:
        if df[column].dtype in ['int64', 'float64']:
            print(f"- {column}")

    print("\nThese columns were text (ignored for statistics):")
    for column in df.columns:
        if df[column].dtype not in ['int64', 'float64']:
            print(f"- {column}")

# Example usage
add_numeric_statistics_to_csv('processed_ecgs_all.csv', 'output_with_stats_7.csv')