import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats


cell_count_filename = '/Users/preethi/Downloads/Teiko_TA/TeikoTechnicalzip/cell-count.csv'
cell_count_df = pd.read_csv(cell_count_filename)
cell_populations = ['b_cell', 'cd8_t_cell', 'cd4_t_cell', 'nk_cell', 'monocyte']

#Question 1
def calculate_relative_frequency(input_df, output_filename):
    # Process each sample
    output = []
    for _, row in input_df.iterrows():
        total = sum(row[cell_populations])
        for population in cell_populations:
            count = row[population]
            percentage = (count / total) * 100 if total > 0 else 0
            output.append({
                'sample': row['sample'],
                'total_count': total,
                'population': population,
                'count': count,
                'percentage': percentage
            })
    
    # Save results
    pd.DataFrame(output).to_csv(output_filename, index=False)

calculate_relative_frequency(cell_count_df, "relative_frequency_output.csv")

#Question 2
#Step 1: Calculating Relative Frequencies
def calculate_relative_frequencies(df):
    # Calculate total cells for each sample
    df['total_count'] = df[cell_populations].sum(axis=1)
    
    # Calculate percentages for each cell population
    for population in cell_populations:
        df[f'{population}_pct'] = (df[population] / df['total_count']) * 100
    
    return df

df = calculate_relative_frequencies(cell_count_df)

# Step 2: Filter for tr1 treatment and PBMC samples only
tr1_df = df[(df['treatment'] == 'tr1') & (df['sample_type'] == 'PBMC')]

# Step 3: Statistical analysis
stats_results = []

for population in cell_populations:
    # Extract data for responders and non-responders
    responders = tr1_df[tr1_df['response'] == 'y'][f'{population}_pct']
    non_responders = tr1_df[tr1_df['response'] == 'n'][f'{population}_pct']
    
    # Perform t-test
    t_stat, p_val = stats.ttest_ind(responders, non_responders)
    
    # Store results
    stats_results.append({
        'population': population.replace('_', ' ').title(),
        't_statistic': round(t_stat, 3),
        'p_value': round(p_val, 4),
        'significant': 'Yes' if p_val < 0.05 else 'No'
    })

stats_df = pd.DataFrame(stats_results)

#Step 4: Generate the boxplots for each population
for population in cell_populations:
    plt.figure(figsize=(10, 7))
    
    # Create plotting dataframe with readable labels
    plot_df = tr1_df.copy()
    plot_df['Group'] = plot_df['response'].map({'y': 'Responder', 'n': 'Non-Responder'})
    
    # Create boxplot with custom styling
    sns.boxplot(x='Group', y=f'{population}_pct', data=plot_df, color='#3EB489')
    
    # Add formatting
    plt.title(f'Relative {population} Cell Count Frequencies in Treatment 1 Responders and Non-Responders')
    plt.ylabel('Relative Frequency (%)')
    plt.grid(axis='y', linestyle='-', alpha=0.2)
    plt.show()
    plt.close()

print("Statistical Analysis Results:")
print(stats_df)

print("\nBased on the Statistical Analysis Results and Boxplots, " \
"Cd4 T cells and Monocytes' relative cell frequencies are significantly " \
"different between responders and non-responders. The boxplots show " \
"Cd4 T cells having higher relative frequencies in responders and " \
"Monocytes having higher relative frequencies in non-responders. " \
"These findings suggest that CD4 T cells and monocytes may play " \
"key roles in predicting response to treatment 1.")