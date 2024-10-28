import pandas as pd


def filter_valid_og():
    # Load the combined data from totalFC.xlsx
    df = pd.read_excel('totalFC.xlsx', header=None)

    # Rename columns to match the required format
    df.columns = ['OGID', 'logFC', 'pvalue', 'genomeID']

    # Convert pvalue and logFC columns to numeric values, coercing errors to NaN
    df['pvalue'] = pd.to_numeric(df['pvalue'], errors='coerce')
    df['logFC'] = pd.to_numeric(df['logFC'], errors='coerce')

    # Filter rows based on the initial conditions for "유효한" genomeID
    valid_conditions = (df['pvalue'] <= 0.05) & (df['logFC'].abs() >= 0.5)
    valid_genomes_df = df[valid_conditions]

    # Group by OGID
    grouped = df.groupby('OGID')

    valid_og = []

    for name, group in grouped:
        total_genomes = group['genomeID'].nunique()
        valid_genomes = valid_genomes_df[valid_genomes_df['OGID'] == name]['genomeID'].nunique()

        if total_genomes >= 4 and valid_genomes >= total_genomes / 2:
            valid_genomes_group = valid_genomes_df[valid_genomes_df['OGID'] == name]
            if (valid_genomes_group['logFC'] > 0).all() or (valid_genomes_group['logFC'] < 0).all():
                valid_og.append(group)

    # Concatenate valid OGs into a single dataframe
    if valid_og:
        valid_df = pd.concat(valid_og, ignore_index=True)
    else:
        valid_df = pd.DataFrame(columns=['OGID', 'logFC', 'pvalue', 'genomeID'])

    # Write the valid dataframe to an Excel file
    valid_df.to_excel('validOG.xlsx', index=False)


# Execute the function
filter_valid_og()
