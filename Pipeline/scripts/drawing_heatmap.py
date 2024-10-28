import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


def transform_data_and_plot_heatmap():
    # Load the combined data from validOG.xlsx
    df = pd.read_excel('validOG.xlsx')

    # Drop the pvalue column
    df = df.drop(columns=['pvalue'])

    # Pivot the dataframe to have genomeID as columns and OGID as index
    df_pivot = df.pivot(index='OGID', columns='genomeID', values='logFC')

    # Save the transformed dataframe to a new Excel file
    df_pivot.to_excel('transformed_validOG.xlsx')

    # Calculate the mean of logFC for each OGID, ignoring NaNs
    df_pivot['mean_logFC'] = df_pivot.mean(axis=1, skipna=True)

    # Sort by the mean_logFC
    df_pivot = df_pivot.sort_values(by='mean_logFC', ascending=False)

    # Drop the mean_logFC column as it's no longer needed
    df_pivot = df_pivot.drop(columns=['mean_logFC'])

    # Define a custom diverging color map
    cmap = sns.diverging_palette(240, 10, n=9, as_cmap=True)  # Blue to red diverging palette

    # Create a heatmap with custom color map
    plt.figure(figsize=(20, 20))  # Adjust the figure size as needed
    sns.heatmap(df_pivot, cmap=cmap, annot=False, linewidths=.5, center=0, cbar_kws={'label': 'logFC'})

    # Set the title and labels
    plt.title('Heatmap of logFC by OGID and genomeID')
    plt.xlabel('genomeID')
    plt.ylabel('OGID')

    # Save the heatmap to a file
    plt.savefig('heatmap_logFC.png')

    # Show the heatmap
    plt.show()


# Execute the function
transform_data_and_plot_heatmap()
